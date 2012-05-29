
import os, collections, datetime, shutil, tarfile

from nesoni import io


def mine_logs(filenames, filter=lambda sample, field: True):
    samples = [ ]
    fields = [ ]
    data = { }
    
    for filename in filenames:
        f = open(filename,'rb')
        for line in f:
            if not line.startswith('(> '): continue
            sample, value, field = line[3:].rstrip().split(None,2)
            if not filter(sample, field):
                continue
            if sample not in samples:
                samples.append(sample)
            if field not in fields:
                fields.append(field)
            data[ (sample,field) ] = value.replace(',','')
        f.close()
    
    records = [ ]
    for field in fields:
        records.append(collections.OrderedDict(
           [ ('Statistic', field) ] + [ (item, data.get((item,field),'')) for item in samples ]
        ))
    return records


STYLE = """
body { font-family: sans-serif; margin-left: 5em; margin-right: 5em; }
img { vertical-align: middle; }
a { text-decoration: none; }
h1,h2,h3,h4,h5,h6 { margin-top: 2em; }
"""

class Reporter(object):
    def __init__(self, directory, title, file_prefix=''):
        self.workspace = io.Workspace(directory, must_exist=False)
        self.file_prefix = file_prefix
        if self.file_prefix: self.file_prefix += '-'
        
        self.f = self.workspace.open('index.html','wb')
    
        print >> self.f, '<html><head>'
        print >> self.f, '<title>%s</title>' % title
        print >> self.f, '<style>%s</style>' % STYLE
        print >> self.f, '</head><body>'
        print >> self.f, '<h1>%s</h1>' % title
        self.p( datetime.date.today().strftime('%e %B %Y') )

    def close(self):
        self.f.close()

    def write(self, text):
        self.f.write(text)

    def heading(self, text):
        print >> self.f, '<h2>%s</h2>' % text

    def p(self, text):
        print >> self.f, '<p>%s</p>' % text

    def href(self, filename, title=None, image=False):
        relative = self.workspace.path_as_relative_path(filename)
        if title is None:
           title = os.path.split(filename)[1]
           
           size = os.stat(filename).st_size
           if size >= 1<<30:
               title += ' (%.1fGb)' % (float(size)/(1<<30))
           elif size >= 1<<20: 
               title += ' (%.1fMb)' % (float(size)/(1<<20))
           elif size >= 1<<10: 
               title += ' (%.1fkb)' % (float(size)/(1<<10))
        
        if image:
           thumb_name = 'thumb-'+relative 
           thumb_filename = self.workspace/thumb_name
           io.execute(['convert', '-thumbnail', '50x50', filename, thumb_filename])
           title = ('<span style="display: inline-block; width: 50px;"><img src="%s"/></span> ' % thumb_name) + title 
           
        return '<a href="%s">%s</a>' % (relative, title)

    def get(self, filename, name=None, title=None, prefix=None, image=False):
        if name is None:
            name = os.path.split(filename)[1]
        if prefix is None:
            prefix = self.file_prefix
        dest = self.workspace / (prefix+name)
        shutil.copyfile(filename, dest)
        return self.href(dest, title, image)
    
    def tar(self, tar_name, filenames, title=None):
        dest = self.workspace / (self.file_prefix+tar_name)
        
        tarf = tarfile.open(dest, 'w:gz')
        for filename in filenames:
            if isinstance(filename, tuple):
                filename, destname = filename
            else:
                destname = os.path.split(filename)[1]            
            tarf.add(filename, destname)        
        tarf.close()
        
        return self.href(dest, title)
        
    def report_logs(self, name, logs, filter=lambda sample, field: True):
        filename = self.workspace / (self.file_prefix + name + '.csv')        
        io.write_csv(filename, mine_logs(logs, filter))        
        self.p(self.href(filename))
        
    def report_heatmap(self, action):
        self.p(
            self.get(action.prefix + '.png', image=True) + ' &sdot; ' +
            self.get(action.prefix + '.csv', title='[spreadsheet]') 
        )




