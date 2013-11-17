
import os, collections, datetime, shutil, tarfile

from nesoni import io, config


def mine_logs(filenames, filter=lambda sample, field: True, commas=False):
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
            if not commas:
                value = value.replace(',','')
            data[ (sample,field) ] = value
        f.close()
    
    records = [ ]
    for field in fields:
        records.append(collections.OrderedDict(
           [ ('', field) ] + [ (item, data.get((item,field),'')) for item in samples ]
        ))
    return records

@config.String_flag('dest')
@config.String_flag('source')
class Copy(config.Action):
    dest = None
    source = None
    
    def ident(self): 
        return super(Copy,self).ident() + '--' + (self.dest or '')
    
    def run(self):
        shutil.copyfile(self.source, self.dest)


@config.Section('files', 'filenames or tuples (filename, destname)')
class Tar(config.Action_with_prefix):
    dest = None
    files = [ ]
    
    def ident(self): 
        return super(Tar,self).ident() + '--' + (self.dest or '')
    
    def run(self):
        tarf = tarfile.open(self.prefix+'.tar.gz', 'w:gz')
        for filename in self.files:
            if isinstance(filename, tuple):
                filename, destname = filename
            else:
                destname = os.path.split(filename)[1]            
            tarf.add(filename, os.path.join(os.path.split(self.prefix)[1],destname))
        tarf.close()


STYLE = """
body { font-family: sans-serif; margin-left: 5em; margin-right: 5em; }
img { vertical-align: middle; }
a { text-decoration: none; }
h1,h2,h3,h4,h5,h6 { margin-top: 2em; }
td,th { padding: 0.5em; }
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

    def subheading(self, text):
        print >> self.f, '<h3>%s</h3>' % text

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
           thumb_name = 'thumb-'+os.path.splitext(relative)[0]+'.png'
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
        
        #Copy(
        #    dest = dest,
        #    source = filename,
        #    ).make()
        io.symbolic_link(source=filename, link_name=dest)
        
        return self.href(dest, title, image)
    
    def tar(self, tar_name, filenames, title=None):
        dest = self.workspace / (self.file_prefix+tar_name)
        
        Tar(
            dest,
            files = filenames,
        ).make()        
        #tarf = tarfile.open(dest, 'w:gz')
        #for filename in filenames:
        #    if isinstance(filename, tuple):
        #        filename, destname = filename
        #    else:
        #        destname = os.path.split(filename)[1]            
        #    tarf.add(filename, destname)        
        #tarf.close()
        
        return self.href(dest+'.tar.gz', title)
        
    def report_logs(self, name, logs, filter=(lambda sample, field: True), renaming={}):
        table = mine_logs(logs, filter)
        
        if name:
            filename = self.workspace / (self.file_prefix + name + '.csv')        
            io.write_csv(filename, table)
            self.p(self.href(filename))
        
        if table:
            table = mine_logs(logs, filter, commas=True)
            
            self.write('<table style="font-size: 50%">\n')
            self.write('<tr>\n')
            for key in table[0].keys():
                self.write('<th>'+key+'</th>')
            self.write('</tr>\n')
            for row in table:
                self.write('<tr>\n')
                for i,value in enumerate(row.values()):
                    if i == 0:
                        value = renaming.get(value,value)
                    self.write('<td>'+value+'</td>')
                self.write('</tr>\n')
            self.write('</table>\n')
                
        
    def report_heatmap(self, action, has_csv=True):
        self.p(
            self.get(action.prefix + '.png', image=True) + 
            (' &sdot; ' +self.get(action.prefix + '.csv', title='[spreadsheet]') if has_csv else '')
        )
    
    #def report_test(self, action):
    #    prefix = os.path.basename(action.prefix)
    #    image = self.get(action.prefix + '-heatmap.png', image=True, title='')
    #    sig_csv = self.get(action.prefix + '.csv', title='[DE genes table]')
    #    all_csv = self.get(action.prefix + '-all.csv', title='[All genes table]')
    #    maybe_maplot = (
    #        (' &sdot; ' + self.get(action.prefix+'.png', image=True,title='[MA-plot]'))
    #        if os.path.exists(action.prefix+'.png') else ''
    #        )
    #    info = self.get(action.prefix + '-info.txt', title='[Info]')        
    #    
    #    text = (
    #        '<table><tr>\n'
    #        '<td valign="top">%(image)s</td>\n'
    #        '<td valign="top"><b>%(prefix)s</b>\n'
    #        '<br/>%(sig_csv)s &sdot; %(all_csv)s &sdot; %(info)s %(maybe_maplot)s\n'
    #        '</tr></table>'
    #        ) % locals()
    #    
    #    self.write(text)
        
        #self.p(
        #    self.get(action.prefix + '-heatmap.png', image=True, title='[heatmap]') +
        #    (' &sdot; ' + self.get(action.prefix+'.png', image=True,title='[MA-plot]')
        #        if os.path.exists(action.prefix+'.png') else '') +
        #    ' &sdot; ' +
        #    self.get(action.prefix + '.csv') +
        #    ' &sdot; ' +
        #    self.get(action.prefix+'-info.txt', title='[info]')
        #)












