
import os, collections, datetime, shutil, zipfile

from nesoni import io


def mine_logs(filenames, omit=[]):
    samples = [ ]
    fields = [ ]
    data = { }
    
    for filename in filenames:
        f = open(filename,'rb')
        for line in f:
            if not line.startswith('(> '): continue
            sample, value, field = line[3:].rstrip().split(None,2)
            if field in omit:
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



class Reporter(object):
    def __init__(self, directory, title, file_prefix=''):
        self.workspace = io.Workspace(directory, must_exist=False)
        self.file_prefix = file_prefix
        if self.file_prefix: self.file_prefix += '-'
        
        self.f = self.workspace.open('index.html','wb')
    
        print >> self.f, '<html><head>'
        print >> self.f, '<title>%s</title>' % title
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

    def href(self, filename, title=None):
        relative = self.workspace.path_as_relative_path(filename)
        if title is None:
           title = os.path.split(filename)[1] 
        return '<a href="%s">%s</a>' % (relative, title)

    def get(self, filename, name=None, title=None, prefix=None):
        if name is None:
            name = os.path.split(filename)[1]
        if prefix is None:
            prefix = self.file_prefix
        dest = self.workspace / (prefix+name)
        shutil.copyfile(filename, dest)
        return self.href(dest, title)
    
    def zip(self, zip_name, filenames, title=None):
        dest = self.workspace / (self.file_prefix+zip_name)
        
        zipf = zipfile.ZipFile(dest, 'w', zipfile.ZIP_DEFLATED, True)
        for filename in filenames:
            if isinstance(filename, tuple):
                filename, destname = filename
            else:
                destname = os.path.split(filename)[1]            
            zipf.write(filename, destname)        
        zipf.close()
        
        return self.href(dest, title)
        
    def report_logs(self, name, logs, omit=[]):
        filename = self.workspace / (self.file_prefix + name + '.csv')        
        io.write_csv(filename, mine_logs(logs, omit))        
        self.p(self.href(filename))
        
    def report_heatmap(self, action):
        self.p(
            self.get(action.prefix + '.png') + ' &nbsp; ' +
            self.get(action.prefix + '.csv', title='[spreadsheet]') 
        )




