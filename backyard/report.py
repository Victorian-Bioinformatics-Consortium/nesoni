
import os, sys, datetime
from os.path import join

from nesoni import grace, io

HEAD = """
<html>
<head>
<title>%(title)s</title>
</head>
<body>
<h1>%(title)s</h1>
<p>%(today)s</p>
"""

TAIL = '</body></html>'

def p(f, text):
    print >> f, '<p>%(text)s</p>' % locals()

def bullet(f, text):
    print >> f, '<ul><li>%(text)s</li></ul>' % locals()

def section(f, title):
    print >> f, '<h2 style="clear: both;">%(title)s</h2>' % locals()

def subsection(f, title):
    print >> f, '<div style="float: left; margin: 0.5em; padding: 0.5em; border: 1px solid black;"><h3>%(title)s</h3>' % locals()

def end_subsection(f):
    print >> f, '</div>'

def pre(f, text):
    print >> f, '<pre>%(text)s</pre>' % locals()

#def copy(src, dest): #TODO: escaping
#    print >> sys.stderr, 'copy '+src+' -> '+dest
#    assert 0 == os.system('cp %(src)s %(dest)s' % locals())

def extract(filename, start_condition):
    f = open(filename,'rU')
    result = [ ]
    started = False
    for line in f:
        if not started and start_condition(line):
            started = True
        if started:
            result.append(line)
    return ''.join(result)

def report_main(args):
    title, args = grace.get_option_value(args, '--title', str, 'Report')
    short_name, args = grace.get_option_value(args, '--short', str, 'files')
    show_refalign, args = grace.get_option_value(args, '--show-refalign', grace.as_bool, True)

    output_dir, args = args[0], args[1:]
    
    reference_filenames = [ ]
    clip_filenames = [ ]
    align_dirs = [ ]
    count_log_filenames = [ ]
    extra_items = [ ]
    extra_files = [ ]
    def file(args): extra_files.append((args[0], ' '.join(args[1:])))
    def extra(args): extra_items.extend(args)
    def reference(args): reference_filenames.extend(args)
    def clips(args): clip_filenames.extend(args)
    def aligns(args): align_dirs.extend(args)
    def count_log(args): count_log_filenames.extend(args)
    
    grace.execute(args, [reference, clips, aligns, extra, file, count_log])
        
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    
    file_dir = join(output_dir, short_name)
    if not os.path.isdir(file_dir): os.mkdir(file_dir)
    for item in os.listdir(file_dir):
        os.unlink(join(file_dir, item))
        
    for filename in reference_filenames:
        io.copy_file(filename, join(file_dir, os.path.basename(filename)))
    
    for filename, desc in extra_files:
        io.copy_file(filename, join(output_dir, os.path.basename(filename)))

    pairs = False
    for directory in align_dirs:
        name = os.path.basename(directory)
        io.copy_file(join(directory,'report.txt'), join(file_dir, name + '-report.txt'))
        for extension in [
            '-depth.userplot',
            '-ambiguous-depth.userplot',
            '-pairspan-depth.userplot',
            '-ambiguous-pairspan-depth.userplot',
        ]:
            filenames = [ item for item in os.listdir(directory)
                          if item.endswith(extension)
                          and not item.endswith('-ambiguous'+extension)
                          and not item.endswith('-pairspan'+extension) ]
            for filename in filenames:
                if len(filenames) == 1:
                    dest = name + extension
                else:
                    dest = name + '-' + filename
                io.copy_file(join(directory,filename), join(file_dir, dest))
                
                if 'pairspan' in extension: pairs = True
    
    today = datetime.date.today().strftime('%e %B %Y')
    
    f = open(join(output_dir, 'index.html'),'wb')
    print >> f, HEAD % locals()
    
    section(f, 'Results')
    
    for item in extra_items:
        p(f, item)
    
    for filename, desc in extra_files:
        name = os.path.basename(filename)
        p(f, '<a href="%(name)s">%(name)s</a> - %(desc)s' % locals())
    
    p(f, '<a href="%(short_name)s.zip">%(short_name)s.zip</a>' % locals())
    
    for filename in reference_filenames:
        bullet(f, os.path.basename(filename) + ' - reference')
    
    bullet(f, '...-report.txt - report on SNPs and indels found')
    
    p(f,'Different kinds of userplot:')
    bullet(f,'...-depth.userplot - depth of coverage of unambiguously aligned reads')
    bullet(f,'...-ambiguous-depth.userplot - depth of coverage, including reads that hit multiple locations')
    if pairs:
        bullet(f,'...-pairspan-depth.userplot - depth, including the space between reads in read-pairs')
        bullet(f,'...-ambiguous-pairspan-depth.userplot - as above, but including reads that hit multiple locations')
    
    if clip_filenames:
        section(f, 'Read clipping')
        
        for filename in clip_filenames:
            assert filename.endswith('_log.txt')
            name = os.path.basename(filename[:-8])
            text = extract(filename, lambda line: line.startswith('Fragments:') or line.startswith('Single reads') or line.startswith('Pairs'))
            subsection(f, name)
            pre(f, text)
            end_subsection(f)            

    if count_log_filenames:
        section(f, 'Counting alignments to genes')
        
        for filename in count_log_filenames:
            pre(f, open(filename,'rb').read())

    if align_dirs and show_refalign:
        section(f, 'Reference alignment')
        for directory in align_dirs:
            name = os.path.basename(directory)
            text = extract(join(directory, 'consensus_log.txt'),
                           lambda line: 'reads/pairs' in line or 'unmapped' in line)
            text = text.replace('(discarded)','')
            text = text.replace('reads/pairs kept', 'aligned unambiguously')
            subsection(f, name)
            pre(f, text)
            end_subsection(f)
    
    print >> f, TAIL % locals()
    f.close()
    
    zip_filename = join(output_dir, short_name + '.zip')
    if os.path.exists(zip_filename): os.unlink(zip_filename)
    assert 0 == os.system('cd %(output_dir)s ; zip %(short_name)s.zip %(short_name)s/* ' % locals())
    for item in os.listdir(file_dir):
        os.unlink(join(file_dir, item))
    os.rmdir(file_dir)
