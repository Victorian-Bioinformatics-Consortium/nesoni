
__all__ = """
   coordinator
   generate 
   remake_needed remake_clear do_nothing
   future parallel_imap parallel_map parallel_for
   Stage process barrier stage stage_function
   make process_make 
   Execute Make
   configure_making run_script run_tool run_toolbox
""".strip().split()

import multiprocessing 
from multiprocessing import managers
import threading, sys, os, signal, atexit, time, base64, socket, warnings
import cPickle as pickle

from nesoni import grace, config

class Child_exception(grace.Error): pass


def chunk(iterable, chunk_size):
    items = [ ]
    for item in iterable:
        items.append(item)
        if len(items) >= chunk_size:
            yield items
            items = [ ]
    if items:
        yield items

def interleave(iterators):
    i = 0
    iterators = list(iterators)
    while iterators:
        i = i%len(iterators)
        try:
            yield iterators[i].next()
            i = i+1
        except StopIteration:
            del iterators[i]




def process_identity():
    return (socket.gethostname(), os.getpid())


# Manager/coordinator ===================


    
class My_coordinator:
    """ LIFO allocation of cores
        - LIFO is generally more memory efficient
        - Behaves as expected when only one core used
        
        Processes are assumed to start owning one core.
        
        Keep track of the number of cores a process has is it's own business.
        
        Maintaining LIFOness makes this is a delicate and subtle dance.
        The implementation of future(...) below is non-obvious.      
    """
    def __init__(self):
        self.lock = threading.RLock()
        self.waiters = [ ]
        self.cores = multiprocessing.cpu_count()
        self.used = 1 #Main process
        
        self.pids = set()
        self.statuses = { }

    def time(self):
        """ A common source of timestamps, if spread over many nodes. """
        return time.time()

    def add_pid(self, pid):
        self.pids.add(pid)
    
    def remove_pid(self, pid):
        self.pids.remove(pid)
    
    def kill_all(self):
        with self.lock:
            if self.pids:
                print >> sys.stderr, 'Killing', ' '.join(map(str,sorted(self.pids)))
                while self.pids:
                    os.kill( self.pids.pop(), signal.SIGTERM )

    def _update(self):
        with self.lock:
            for i in xrange(len(self.waiters)-1,-1,-1):
                if self.waiters[i][0] + self.used <= self.cores:
                    self.used += self.waiters[i][0]
                    self.waiters[i][1].set()
                    del self.waiters[i]
    
    def set_cores(self, n):
        with self.lock:
            self.cores = n
            self._update()

    def get_cores(self):
        with self.lock:
            return self.cores
    
    def change_cores_used(self, delta):
        """ Advise of use of more or less cores, no delay allowed. """
        with self.lock:
            self.used += delta
            self._update()
    
    def release_core(self):
        self.change_cores_used(-1)

    def trade_cores(self, old, new):
        with self.lock:
            assert new <= self.cores, 'Don\'t have %d cores.' % new
            self.used -= old
            if self.used + new <= self.cores:
                self.used += new
                return
            event = threading.Event()
            self.waiters.append((new, event))
            self._update()
        event.wait()

    def acquire_core(self):
        self.trade_cores(0,1)

    def set_status(self, identity, value):
        if value:
            self.statuses[identity] = value
        elif identity in self.statuses:
            del self.statuses[identity]
        
        if sys.stderr.isatty():
            items = [ self.statuses[item] for item in sorted(self.statuses) ]
            alloc = 200 / max(1,len(items))
            status = ''
            for item in items:
                if len(item) > alloc:
                    status += item[:alloc]+'> '
                else:
                    status += item + ' '

            #Show in terminal title
            sys.stderr.write('\x1b]2;'+status+'\x07')
            sys.stderr.flush()
            
        

# The coordinator in the manager-process
_COORDINATOR = None
def _get_coordinator():
    global _COORDINATOR
    if not _COORDINATOR:
        _COORDINATOR = My_coordinator()
    return _COORDINATOR


class My_manager(managers.SyncManager):
    class _Server(managers.SyncManager._Server):
        def serve_forever(self):
            signal.signal(signal.SIGINT, signal.SIG_IGN)
            managers.SyncManager._Server.serve_forever(self)
    

My_manager.register('get_coordinator', callable=_get_coordinator)

_MANAGER = None
_AUTHKEY = None
def manager(address=('127.0.0.1',0),authkey=None,connect=False):
    global _MANAGER, _AUTHKEY
    if _MANAGER is None:
        if authkey is None:
            with open('/dev/urandom','rb') as f:
                authkey = base64.b16encode(f.read(256))
        _MANAGER = My_manager(address=address, authkey=authkey)
        _AUTHKEY = authkey
        if connect:
            _MANAGER.connect()
        else:
            _MANAGER.start()            
            atexit.register( coordinator().kill_all )
    return _MANAGER


# Local proxy of the coordinator in the manager-process
_COORDINATOR_PROXY = None
def coordinator():
    global _COORDINATOR_PROXY
    if _COORDINATOR_PROXY is None:
        _COORDINATOR_PROXY = manager().get_coordinator()
    return _COORDINATOR_PROXY


def _in_process(address, authkey, func, args, kwargs):
    global _MANAGER, _COORDINATOR_PROXY
    # Might have weirdness from forking
    _MANAGER = None
    _COORDINATOR_PROXY = None    
    manager(address,authkey,True)
    coordinator().add_pid( os.getpid() )
    try:
        return func(*args,**kwargs)
    finally:
        coordinator().remove_pid( os.getpid() )
        
def start_process(func, *args, **kwargs):
    man = manager()
    p = multiprocessing.Process(target=_in_process,args=(man.address,_AUTHKEY,func,args,kwargs))
    p.start()
    return p

def subprocess_Popen(*args, **kwargs):
    # Force manager start before spawning any subprocesses
    # so it doesn't inhert any pipes
    manager()
    
    import subprocess
    return subprocess.Popen(*args, **kwargs)



# =======================================

class Stage(object):
    """ All the world's a stage,
        And all the men and women merely players 
    
        Use this class to synchronize with sets of 
        processes that you start.
        
        Example:
        
        
        stage = Stage()
        stage.process(my_func1,...)
        stage.process(my_func2,...)
        ...
        
        stage.barrier()   #Wait for all processes started by this stage to finish

        
        Limitations: 
        A stage object can not be passed to a different process.
    """
    def __init__(self):
        self.futures = [ ]

    def add(self, future):
        """ Add an existing process to this stage's collection. 
        
            Actually, it can be anything that can be called with
            no arguments.
        """
        if not self.futures:
            LOCAL.stages.add(self)
        self.futures.append(future)
        
    def process(self, func, *args, **kwargs):
        """ Create a new process that will execute func, and
            add it to this stage's collection. """
        item = future(func, *args, **kwargs)
        self.add(item)
        return item

    def barrier(self):
        """ Wait for all processes that have been added to this stage
            to finish. """
        if self.futures:
            LOCAL.stages.remove(self)

            exceptions = [ ]
            while self.futures:
                try:
                    self.futures.pop()()
                except Exception, e:
                    exceptions.append(e)
            if exceptions:
                raise Child_exception('Child failures: %d' % len(exceptions), exceptions)





def set_locals():
    global LOCAL
    LOCAL = threading.local()
    LOCAL.abort_make = False
    LOCAL.do_nothing = False 
    LOCAL.time = 0 #Note: do not run this code before the year 1970
    LOCAL.stages = set() #Stages with processes in them, so we can warn if they don't have .barrier() called on them
    LOCAL.stage = Stage() #Default stage. Deprecated.
set_locals()

def _check_stages():
    if LOCAL.stages:
        warnings.warn('Exited without calling .barrier() on all Stages.')  
atexit.register(_check_stages)


def remake_needed():
    """ Force all tools to be re-run.     
    """
    LOCAL.time = coordinator().time()

def remake_clear():
    """ Indicate that any following tool invocation do not depend
        on previous tool invocations. """
    LOCAL.remake = 0

def abort_makes():
    """ Immediately abort if any tool needs to be run.
    
        Allows checking of what would be run without actually running it.    
    """
    LOCAL.abort_make = True

def do_nothing():
    LOCAL.do_nothing = True




def _run_future(time,abort_make,do_nothing, func, args, kwargs, sender):
    set_locals()
    LOCAL.time = time
    LOCAL.abort_make = abort_make
    LOCAL.do_nothing = do_nothing
    result = None
    exception = None
    try:
        result = func(*args, **kwargs)
        barrier()        
        assert not LOCAL.stages, 'Process completed without calling .barrier() on all Stages.'
    except object, e:
        config.report_exception()
        exception = e
    
    #Send result to parent before releasing core
    sender.send((LOCAL.time, exception, result))

    #Give core back to parent
    coordinator().release_core()

    sender.close()


def future(func, *args, **kwargs):
    """
    Underlying synchronization mechanism.
    
    Create a new process to run a function.
    
    The return value can be later called with
    no arguments to get the result of the function. 
    This has the side effect of synchronizing with 
    the process that was created, and this process
    being "infected" with the need to remake if 
    necessary.
    
    Limitation:
    A future can not be passed to a different process.
    """

    # Pause, let your mind expand.
    # What follows is as it must be,
    # unless I have made a mistake.    
    
    receiver, sender = multiprocessing.Pipe(False)

    #Give core to process we start
    p = start_process(_run_future,LOCAL.time,LOCAL.abort_make,LOCAL.do_nothing,func,args,kwargs,sender)
    
    #Get another for ourselves
    coordinator().acquire_core()

    received = [ ]
    def local_func():
        if not received:
            #Get the core back from the process we started
            if not receiver.poll():
                # Don't need our core while we're waiting
                coordinator().release_core()
                receiver.poll(None)
                coordinator().acquire_core()
            received.append( receiver.recv() )            
            receiver.close()
            p.join()
        
        time, exception, result = received[0]        
        LOCAL.time = max(LOCAL.time, time)
        if exception:
            raise Child_exception('Child failed', exception)
        return result
    return local_func


def parallel_imap(func,iterable):
    # This may be memory inefficient for long iterators
    return (item() for item in [ future(func,item2) for item2 in iterable ])

def parallel_map(func, iterable):
    return list(parallel_imap(func, iterable))

def parallel_for(iterable):
    """ Execute a "for loop" in parallel.
        Use this as a function decorator. 
    """
    def thread_for(func):
        for item in parallel_imap(func, iterable):
            pass
    return thread_for


def process(func, *args, **kwargs):
    """ Deprecated. Use Stage objects.
    
        Start a new process. 
    """
    return LOCAL.stage.process(func, *args, **kwargs)

#def thread(func, *args, **kwargs):
#    LOCAL.parallels.append(future(func, *args, local=True, **kwargs))
#    return func

def barrier():
    """ Deprecated. Use Stage objects.
    
        Wait for all processes started by this process to finish.    
        (Except for any processes explicitly put in a stage.)
    """
    LOCAL.stage.barrier()


def stage(func, *args, **kwargs):
    """ Deprecated. Use Stage objects.
    
        Call a function, and wait for all processes started by it to finish    
        Can be used as a function decorator, but probably shouldn't be.
    """
    old = LOCAL.stage
    LOCAL.stage = Stage()
    result = None
    try:
        result = func(*args, **kwargs)
    finally:
        barrier()    
        LOCAL.stage = old
    return result


def stage_function(func):
    """ Deprecated. Use Stage objects.
    
        Ensure processes started by a function or method complete
        before the function returns.
    """
    def inner(*args, **kwargs):
        return stage(func, *args, **kwargs)
    return inner


    




# Make ===============================================

def _get_timestamp(action):
    """ Look for ident() in .state subdirectory of current directory.
        If pickled value matches return the timestamp.
    """
    try:
        if not os.path.exists('.state'):
            os.mkdir('.state')
    
        filename = os.path.join('.state', grace.filesystem_friendly_name(action.ident()))
        if os.path.exists(filename):
            with open(filename,'rb') as f:
                old = pickle.load(f)
            
            if action == old:
                if not hasattr(old, 'timestamp'):
                    return None                        
                return old.timestamp
            
            #for parameter in self.parameters:
            #    if parameter.get(self) != parameter.get(old):
            #        print >> sys.stderr, parameter.name, parameter.get(old), '->', parameter.get(self)
            
    except Exception, error:
        import traceback
        traceback.print_exc()
        print >> sys.stderr, 'Error making %s, re-running: %s' % (action.ident(), error)

    return None
    
def _run_and_save_state(action, timestamp):
    filename = os.path.join('.state', grace.filesystem_friendly_name(action.ident()))
    temp_filename = os.path.join('.state', 'temp-' + grace.filesystem_friendly_name(action.ident()))
    
    if os.path.exists(filename):
        os.unlink(filename)
    
    if LOCAL.do_nothing:
        result = None
    else:
        result = action.run()
    
    LOCAL.time = max(LOCAL.time, timestamp)
    action.timestamp = timestamp
    with open(temp_filename,'wb') as f:
        pickle.dump(action, f)
    os.rename(temp_filename, filename)
    
    return result

def _make_inner(action):
    timestamp = coordinator().time()
    assert timestamp > LOCAL.time, 'Time running in reverse.'
    
    cores = action.cores_required()
    if cores > 1:
        coordinator().trade_cores(1, cores)
    try:        
        config.write_colored_text(sys.stderr, '\n'+action.describe()+'\n')
        
        if LOCAL.abort_make:
            raise grace.Error('%s would be run. Stopping here.' % action.ident())
        
        grace.status(action.ident())
        try:
            _run_and_save_state(action, timestamp)
        finally:
            grace.status('')
    finally:
        if cores > 1:
            coordinator().trade_cores(cores, 1)
    


def make(action):
    """ Run a tool (an instance of a subclass of nesoni.config.Action) if necessary.
    """
    timestamp = _get_timestamp(action)    
    if timestamp is not None and timestamp >= LOCAL.time:
        LOCAL.time = timestamp
    else:
        _make_inner(action)


def _time_advancer(timestamp):
    def time_advancer():
        LOCAL.time = max(LOCAL.time, timestamp)
    return time_advancer

def process_make(action, stage=None):
    """ This is just a more efficient version of stage.process(make, <action>)
    """
    if stage is None:
        stage = LOCAL.stage
    timestamp = _get_timestamp(action)    
    if timestamp is not None and timestamp >= LOCAL.time:
        stage.add(_time_advancer(timestamp))
    else:
        stage.process(_make_inner,action)



def generate(func, *args, **kwargs):
    """ Run an iterator in a separate process.
    
        For example:
        
          for item in thing_maker(param):
              ...
            
        could be rewritten:
        
          for item in generate(thing_maker, param):
              ...    
    """
    #return iter(func(*args,**kwargs))
    
    receiver, sender = multiprocessing.Pipe(False)
    
    def sender_func():
        try:
            for items in chunk(func(*args,**kwargs), 1<<8):
                sender.send(items)
        except:
            sender.send('error')
            import traceback
            print >> sys.stderr, traceback.format_exc()            
        else:    
            sender.send(None)
        finally:
            coordinator().change_cores_used(-1)
        sender.close()
    
    coordinator().change_cores_used(1)
    process = start_process(sender_func)
    
    def generator():
        while True:
            data = receiver.recv()
            if data is None: break
            if data == 'error': raise Child_exception()
            for item in data:
                yield item
        process.join()
        receiver.close()
    
    return generator()




@config.help("""\
Execute a shell command, optionally reading stdin from a file, \
and optionally sending stdout to a file.
""")
@config.Int_flag('cores','Advise how many cores the command will use.')
@config.Main_section('command','Command to execute', allow_flags=True, empty_is_ok=False)
class Execute(config.Action_filter):
    cores = 1
    command = [ ]

    def cores_required(self):
        return self.cores
    
    def run(self):
        from nesoni import io
        
        f_in = self.begin_input()
        f_out = self.begin_output()
        try:
            io.execute(self.command, stdin=f_in, stdout=f_out)
        finally:
            self.end_output(f_out)
            self.end_input(f_in)



@config.Int_flag('make_cores', 'Approximate number of cores to use.')
@config.Bool_flag('make_force', 'Force everything to be recomputed.')
@config.Bool_flag('make_show', 'Show the first actions that would be made, then abort.')
@config.Bool_flag('make_done', 
    'Do nothing, but mark all actions as done. '
    'This might be useful if there is a trivial parameter change you don\'t want to re-run. '
    'To re-run from a particular point, use this option then delete files from .state/ as needed.')
#@config.String_flag('manager_address', 'IP address of the network interface you want the job manager to listen to.')
#@config.String_flag('launch', 'Command to launch a new python.')
class Make(config.Action):
    make_cores = multiprocessing.cpu_count()
    make_force = False
    make_show = False
    make_done = False
    
#    manager_address = '127.0.0.1'
#    launch = sys.executable
    
    def _before_run(self):
        coordinator().set_cores(self.make_cores)
        if self.make_force:
            remake_needed()
        if self.make_show:
            abort_makes()
        if self.make_done:
            do_nothing()

    def run(self):
        pass


@config.help("""\
Execute a script.
""")
@config.Hidden('function', 'Function to execute.')
@config.Main_section('script_parameters', 'Script parameters.', allow_flags=True)
class Make_script(Make):
    function = None
    script_parameters = [ ]

    def run(self):
        stage(self.function, *self.script_parameters)


def configure_making(args):
    """ Configure make options, return remaining arguments.
    
        Exits program if options are invalid.
    """
    
    try:
        maker = Make()
        leftovers = maker.parse_partial(args)
        if leftovers != args:
            config.write_colored_text(sys.stderr, '\n'+maker.describe('Make options')+'\n')        
        maker.run()    
        return leftovers
    except:
        config.report_exception()
        sys.exit(1)
    

def run_script(function):
    """ Run a workflow script. Various command line flags are parsed,
        then any remaining command line parameters are passed to the 
        function. 
        
        Intended usage:

    
        from nesoni import *
    
        def my_script():
            ...
        
        if __name__ == '__main__':
            run_script(my_script)
    """
    maker = Make_script(function=function)
    config.shell_run(maker, sys.argv[1:], sys.executable + ' ' + sys.argv[0])


def run_tool(action_class):
    """
    Provide a command line interface for an Action.
    """
    args = configure_making(sys.argv[1:])
    config.shell_run(action_class(), args, sys.argv[0])


def run_toolbox(action_classes, script_name=''):
    """
    Provide a command line interface for a list of Actions.
    
    Note:    
    strings included in the action_classes list will be printed 
    as help text, for example to display section headings.
    """
    args = configure_making(sys.argv[1:])
    
    commands = { }
    help = [ '\n' ]
    for item in action_classes:
        if isinstance(item, str):
            help.append(config.wrap(item, 70) + '\n\n')
            continue
        name = item.shell_name()
        commands[ name ] = item
        help.append('    %s\n' % config.colored(1,name+':'))
        help.append(config.wrap(item.help_short, 70, '        ') + '\n\n')

    help.append('\nMake options:\n'+Make().describe('', show_help=True)+'\n')

    if not args:
        config.write_colored_text(sys.stdout, ''.join(help)+'\n\n')
        sys.exit(1)
        
    try:        
        command, args = args[0], args[1:]
        
        mangled_command = command.lower().rstrip(':')
        if mangled_command not in commands:
            raise grace.Error("Don't know how to "+command)        
    except:
        config.report_exception()
        sys.exit(1)

    config.shell_run(commands[mangled_command](), args, (script_name+' ' if script_name else '') + mangled_command+':')


