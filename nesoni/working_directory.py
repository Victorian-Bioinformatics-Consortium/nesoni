
import os, re

from nesoni import grace, io, config, reference_directory

def _is_bool(item):
    return type(item) is bool

def matches(expression, tags):
    priority = { '[':3, ']':0, ':':2, '+':1, '-':1 }
    exp = re.split('(['+''.join('\\'+item for item in priority)+'])',expression)
    exp = [ item for item in exp if item ]
    exp = [ item in tags if item not in priority else item for item in exp ]
    i = 0
    while i < len(exp)-2:
        print exp
        if i == len(exp)-3:
           p = 0
        else:
           p = priority.get(exp[i+3],0)
           
        if exp[i] == '[' and _is_bool(exp[i+1]) and exp[i+2] == ']':
            exp[i:i+3] = [ exp[i+1] ]
            i = 0
        elif p <= 2 and _is_bool(exp[i]) and exp[i+1] == ':' and _is_bool(exp[i+2]):
            exp[i:i+3] = [ exp[i] and exp[i+2] ]
            i = 0
        elif p <= 1 and _is_bool(exp[i]) and exp[i+1] == '+' and _is_bool(exp[i+2]):
            exp[i:i+3] = [ exp[i] or exp[i+2] ]
            i = 0
        elif p <= 1 and _is_bool(exp[i]) and exp[i+1] == '-' and _is_bool(exp[i+2]):
            exp[i:i+3] = [ exp[i] and not exp[i+2] ]
            i = 0
        else:
            i += 1
    assert len(exp) == 1, 'Could not parse "%s"' % expression
    return exp[0]


class Working(io.Workspace):
    def set_reference(self, path):
        self.update_param(reference=self.path_as_relative_path(path))

    def setup_reference(self, filenames, bowtie=False):
        if len(filenames) == 1 and os.path.isdir(filenames[0]):
            self.set_reference(filenames[0])
            return
        
        path = self / 'reference'
        reference_directory.Make_reference(path, filenames=filenames, bowtie=bowtie).run()
        self.set_reference(path)

    def get_reference(self):
        if 'reference' in self.param:
            path = self.relative_path_as_path(self.param['reference'])
        else:
            path = self.working_dir
        
        return reference_directory.Reference(path, must_exist=True)

    def get_filtered_sorted_bam(self):
        filename = self / 'alignments_filtered_sorted.bam'
        assert os.path.exists(filename), 'Alignments in %s haven\'t been filtered, need to run "nesoni filter:" or "nesoni consensus:".' % self.name
        return filename

    def get_tags(self):
        return [self.name,'all']+self.param.get('tags',[])

    def matches(self, expression):
        return matches(expression, self.get_tags())



@config.help("""\
Label a working directory with a list of tags.

(The list is stored in the file <workding_dir>/parameters.)
""")
@config.Main_section('tags', 'Tags to give working directory.\n(Any existing tags are discarded.)')
class Tag(config.Action_with_working_dir):
    tags = [ ]
    
    _workspace_class = Working
    
    def run(self):
        for tag in tags:
            for char in '+-:, \t\'\"':
                assert char not in tag, 'Tags shouldn\'t contain "'+char+'".'
        
        workspace = self.get_workspace()
        workspace.update_param(tags=self.tags)




