
import os, re

from nesoni import grace, io, config, reference_directory

def _is_bool(item):
    return type(item) is bool


def matches(expression, tags):
    tokens = list('[]-:/^')
    exp = re.split('(['+''.join('\\'+item for item in tokens)+'])',expression)
    exp = [ item for item in exp if item ]
    
    def parse2(exp):
        assert exp
        if exp[0] == '[':
           value, exp = parse(exp[1:])
           assert exp and exp[0] == ']', 'expected a closing ]'
           return value, exp[1:]
        assert exp[0] not in tokens, 'didn\'t expect '+exp[0]
        return exp[0] in tags, exp[1:]
    
    def parse1(exp):
        assert exp, 'unexpected end of expression'
        if exp[0] == '-':
            value, exp = parse2(exp[1:])
            return not value, exp
        else:
            value, exp = parse2(exp)
            return value, exp
    
    def parse(exp):
        value, exp = parse1(exp)
        while exp and exp[0] in [':','/','^']:
            operator, exp = exp[0], exp[1:]
            value2, exp = parse1(exp)
            if operator == ':':
                value = value and value2
            elif operator == '/':
                value = value or value2
            else:
                value = (not value2 and value) or (not value and value2)
        return value, exp
    
    try:
        value, exp = parse(exp)
        assert not exp, 'don\'t know what to do with: '+''.join(exp)
    except AssertionError, e:
        raise grace.Error('Could not parse: '+expression+', '+e.args[0])
    return value


def select_and_sort(select_expression, sort_expression, items, get_tags=lambda item: item.get_tags()):
    """ Select items based on select_expression then sort by sort_expression. """
    items = [ item for item in items
              if matches(select_expression, get_tags(item)) ]

    if not sort_expression:
        parts = []
    else:
        parts = sort_expression.split(',')
    
    def key(item):
        tags = get_tags(item)
        return [ 0 if matches(part, tags) else 1
                 for part in parts ]

    items.sort(key=key)
    return items


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
        for tag in self.tags:
            for char in '+-:, \t\'\"':
                assert char not in tag, 'Tags shouldn\'t contain "'+char+'".'
        
        workspace = self.get_workspace()
        workspace.update_param(tags=self.tags)




