
from nesoni import grace

import numpy
from numpy import random, linalg

import gtk

from matplotlib.figure import Figure

# uncomment to select /GTK/GTKAgg/GTKCairo
#from matplotlib.backends.backend_gtk import FigureCanvasGTK as FigureCanvas
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
#from matplotlib.backends.backend_gtkcairo import FigureCanvasGTKCairo as FigureCanvas

# or NavigationToolbar for classic
#from matplotlib.backends.backend_gtk import NavigationToolbar2GTK as NavigationToolbar
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar



from matplotlib import artist, transforms, path, patches, font_manager

class Arrow(artist.Artist):
    zorder = 3

    def __init__(self, points, text, up_side = True, **kwargs):
        #self.start = numpy.asarray(start,'float64')
        #self.end = numpy.asarray(end,'float64')
        self.points = numpy.asarray(points, 'float64')
        self.text = text
        self.prop = font_manager.FontProperties(size=15)
        self.up_side = up_side
        artist.Artist.__init__(self, **kwargs)

    def draw(self, renderer):
        gc = renderer.new_gc()
        
        gc.set_linewidth(4)
        gc.set_foreground((0.75,0.0,0.75)) #(0.45,0.7,0.7))
        gc.set_antialiased(True)
        self._set_gc_clip(gc)
        gc.set_alpha(1.0)
        
        transform = self.get_transform()
        
        tpoints = numpy.array([ transform.transform_point(item) for item in self.points ])
        
        if len(tpoints) == 1:
            dirs = numpy.array([[ 0.0, 0.0 ]])
        else:
            dirs = tpoints[1:] - tpoints[:-1]
            dirs = numpy.concatenate([
                [ dirs[0] ],
                dirs[:-1] + dirs[1:],
                [ dirs[-1] ],
            ])
            
        #Smooth a bit
        #for i in xrange(3):
        #    newdirs = dirs.copy()
        #    newdirs[1:] += dirs[:-1]
        #    newdirs[:-1] += dirs[1:]
        #    dirs = newdirs

        lengths = numpy.sqrt( dirs[:,0]*dirs[:,0] + dirs[:,1]*dirs[:,1] )
        zero_lengths = (lengths == 0.0)
        lengths[zero_lengths] = 1.0
        dirs[zero_lengths] = [ 1.0, 0.0 ]
        dirs /= lengths[:,None]
        normals = numpy.transpose([
            -dirs[:,1],
            dirs[:,0]
        ])             
        
        offset = 5        
        
        if len(tpoints) > 1:
            line_points = tpoints + normals * offset
            #mypath = path.Path(line_points)
            
            path_points = [ line_points[0], line_points[0] ]
            for i in xrange(1,len(line_points)):
                 path_points.append( (path_points[-1]+line_points[i])*0.5 )
                 path_points.append( line_points[i] )
            path_points.append( line_points[-1] )
            mypath = path.Path(path_points,[path.Path.MOVETO]+[path.Path.CURVE3]*(len(path_points)-1))
                
                    
            renderer.draw_path(gc, mypath, transforms.IdentityTransform())
            
            arrow_size = offset * 1.6 #0.8
            
            mypath = path.Path(numpy.array([
                 line_points[-1] + (-dirs[-1]+normals[-1])*arrow_size,
                 line_points[-1],   
                 line_points[-1] + (-dirs[-1]-normals[-1])*arrow_size,
            ]))
            renderer.draw_path(gc, mypath, transforms.IdentityTransform())
        
        #start = transform.transform_point(self.start)
        #end = transform.transform_point(self.end)
        #
        #if start[0] != start[0]:
        #    #print 'Couldn\'t draw', self.text, self.start, self.end
        #    return
        #
        #vec = end - start
        #length = numpy.sqrt(numpy.sum(vec*vec))
        #if length > 0.0:
        #    unit_vec = vec / length
        #else:
        #    unit_vec = numpy.array([1.0,0.0])
        #    
        #unit_norm = numpy.array([ -unit_vec[1], unit_vec[0] ])
        #        
        #width = 5.0
        #head_width = width*2.0
        #
        #offset = width * 2.0
        #
        #start += unit_norm * (offset if self.up_side else -offset)
        #end   += unit_norm * (offset if self.up_side else -offset)
        #
        #mid = start + unit_vec * max(0.0, length - head_width)
        #
        #mypath = path.Path(numpy.array([
        #    start + unit_norm * width,
        #    mid + unit_norm * width,
        #    mid + unit_norm * head_width,
        #    end,
        #    mid - unit_norm * head_width,
        #    mid - unit_norm * width,
        #    start - unit_norm * width, 
        #    start + unit_norm * width,
        #]))
        #        
        #renderer.draw_path(gc, mypath, transforms.IdentityTransform(), (0.45,0.7,0.7))

        w,h,d = renderer.get_text_width_height_descent(self.text,self.prop,False)
                
        #if h < length:
        if True:
        #if len(tpoints) == 1 or numpy.var(tpoints[:,0])+numpy.var(tpoints[:,1]) > 20**2: #Hmm
            mid = len(tpoints) // 2
            xy =  tpoints[mid] + normals[mid] * (offset*2) + \
                  dirs[mid] * (w*-0.5)
                  
            
            text_orient = ( numpy.arctan2(dirs[mid,1],dirs[mid,0]) /numpy.pi*180.0 ) % 360.0
            
            if not self.up_side:
                r = text_orient * numpy.pi / 180.0
                xy[0] -= -numpy.sin(r) * h
                xy[1] -=  numpy.cos(r) * h

            if text_orient > 90.0 and text_orient < 270.0:
                r = text_orient * numpy.pi / 180.0
                xy[0] += numpy.cos(r) * w
                xy[1] += numpy.sin(r) * w
                xy[0] += -numpy.sin(r) * h
                xy[1] +=  numpy.cos(r) * h
                text_orient -= 180.0

            
            #wtf
            if renderer.flipy():
                canvasw, canvash = renderer.get_canvas_width_height()
                xy[1] = canvash - xy[1]
            
            gc.set_foreground((1.0,1.0,1.0))
            for ox,oy in [(-1,0),(1,0),(0,-1),(0,1)]:
                renderer.draw_text(gc, xy[0]+ox,xy[1]+oy, self.text, self.prop, text_orient, False)
            
            gc.set_foreground((0.0,0.0,0.0))            
            renderer.draw_text(gc, xy[0],xy[1], self.text, self.prop, text_orient, False)


class Quadleaf(object):
    __slots__ = ('ys','xs','values','total')
    
    def visit(self, top, left, bottom, right, should_descend, visitor):
        for i in xrange(len(self.ys)):
            visitor(self.ys[i],self.xs[i],self.values[i])
    

class Quadnode(object):
    __slots__ = ('y', 'x', 'size', 'total', 'children', 'built') 

    def build(self):
        self.children[0] = make_quadtree(*self.children[0])
        self.children[1] = make_quadtree(*self.children[1])
        self.children[2] = make_quadtree(*self.children[2])
        self.children[3] = make_quadtree(*self.children[3])
        self.built = True

    def visit(self, top, left, bottom, right, should_descend, visitor):
        if top > self.y+self.size or \
           left > self.x+self.size or \
           bottom < self.y-self.size or \
           right < self.x-self.size:
            return
        
        if should_descend(self.size, self.total):
            if not self.built: self.build()
            for child in self.children:
                child.visit(top,left,bottom,right,should_descend,visitor)
        else:
            visitor(self.y, self.x, self.total)
        
    
    def child_id(self, y,x):
        if y < self.y:
            if x < self.x:
                return 0
            else:
                return 1
        else:
            if x < self.x:
                return 2
            else:
                return 3 

def make_quadtree(ys,xs,values):
    if len(ys) < 5:
        result = Quadleaf()
        result.ys = ys
        result.xs = xs
        result.values = values
        result.total = numpy.sum(values,0)
        return result
    
    top = numpy.minimum.reduce(ys)
    bottom = numpy.maximum.reduce(ys)
    left = numpy.minimum.reduce(xs)
    right = numpy.maximum.reduce(xs)
    
    result = Quadnode()
    result.y = numpy.average(ys)
    result.x = numpy.average(xs)
    result.size = max(result.y-top,bottom-result.y,result.x-left,right-result.x)
        
    tl = (ys < result.y) & (xs < result.x)
    tr = (ys < result.y) & (xs >= result.x)
    bl = (ys >= result.y) & (xs < result.x)
    br = (ys >= result.y) & (xs >= result.x)
    
    result.built = False
    result.children = [
        (ys[tl],xs[tl],values[tl]),
        (ys[tr],xs[tr],values[tr]),
        (ys[bl],xs[bl],values[bl]),
        (ys[br],xs[br],values[br]),
    ]
    
    #result.total = numpy.sum([ child.total for child in result.children ], 0)
    result.total = numpy.sum(values, 0)
    return result
    
class Dots(artist.Artist):
    """ Quickly draw a collection of dots.
    """
    
    def __init__(self, ys, xs, colors, sizes, zorder=0, **kwargs):
        n = len(ys)
        values = numpy.concatenate(
           [ colors * sizes[:,None],
             sizes[:,None],
             numpy.ones(n)[:,None] ],
           axis=1)
        
        self.quadtree = make_quadtree(ys, xs, values)
        
        self.zorder = zorder

        artist.Artist.__init__(self, **kwargs)

    def draw(self, renderer):
        gc = renderer.new_gc()
        
        gc.set_linewidth(0)
        gc.set_antialiased(True)
        self._set_gc_clip(gc)
        gc.set_alpha(1.0)
                
        transform = self.get_transform()
        
        #Assumes a simple transform
        a = transform.transform_point((0.0,0.0))
        b = transform.transform_point((1.0,1.0))
        
        unit_size = min( 1.0/abs(b[0]-a[0]), 1.0/abs(b[1]-a[1]) ) #* 0.5
        
        clip_box = gc.get_clip_rectangle().inverse_transformed(transform)

        marker_path = path.Path.unit_regular_polygon(5)
        
        def should_descend(size, value):
            return size >= numpy.sqrt(value[3]/value[4]) * unit_size *0.5
        def visitor(y,x,value):
            color = value[:3] / value[3]
            radius = numpy.sqrt( value[3] / value[4] )
            
            tx, ty = transform.transform_point((x,y))
            
            #print tx, ty, radius

            #mypath = path.Path(numpy.array([
            #    [ tx-radius, ty-radius ],
            #    [ tx+radius, ty-radius ],
            #    [ tx+radius, ty+radius ],
            #    [ tx-radius, ty+radius ],
            #    [ tx-radius, ty-radius ],                
            #]))
            
            rtrans = transforms.BboxTransformTo(transforms.Bbox.from_bounds(
                tx,ty,radius,radius
            ))
            #tpath = rtrans.transform_path(upath)
                
            renderer.draw_path(gc, marker_path, rtrans, tuple(color))
        
        self.quadtree.visit(clip_box.y0,clip_box.x0,clip_box.y1,clip_box.x1, should_descend, visitor)







class Graph:
    def __init__(self, node_names, weights):
        self.names = [ ]
        self.weights = weights
        self.name_to_ident = { }
        self.links = [ ]
        for name in node_names:
            self.name_to_ident[name] = len(self.names)
            self.names.append(name)
            self.links.append([])
        
        self.positions = random.random((len(self.links),2))
        
        self.update_amount = 1.5

    def has(self, name):
        return name in self.name_to_ident

    def link(self, name1, name2, sign, step):
        # sign : 1 = same dir, -1 = change dir
        # step : travel in frame of name1
        ident1 = self.name_to_ident[name1]
        ident2 = self.name_to_ident[name2]
        self.links[ident1].append( (ident2, sign, step) )
        self.links[ident2].append( (ident1, sign, sign*-step) )

    def _improve_layout(self, root, max_n, amount,  idents,signs,distances,travels):
        links = self.links
        end = min(len(links), int(max_n)+1)

        idents[0] = root
        signs[0] = 1
        distances[0] = 0
        travels[0] = 0
        j = 1
        
        queued = set()
        queued.add(root)
        
        i = 0
        while i < j and j < end:
            #if i >= len(idents):
            #    while True:
            #        j = random.randint(len(self.links))
            #        if j not in queued: break                  
            #    queued.add(j)
            #    idents.append(j)
            #    distances.append(distances[-1])
        
            ident = idents[i]
            distance_plus_1 = distances[i] + 1
        
            for ident2, sign2, travel2 in links[ident]:
                if ident2 in queued: continue
                
                queued.add(ident2)
                idents[j] = ident2
                signs[j] = signs[i] * sign2
                distances[j] = distance_plus_1
                travels[j] = travels[i] + travel2*signs[i]
                
                j += 1

            i += 1

        if i <= 1: return
        
        end = min(end, j)
        idents = idents[:end]
        signs = signs[:end]
        distances = distances[:end]
        travels = travels[:end]

        node_weights = numpy.array([ self.weights[i] for i in idents ], dtype='float64')
        
        this_positions = self.positions[idents]
                
        offsets = this_positions - self.positions[root][None,:]
        
        #print travels 
        #mainline = (numpy.absolute(travels) == distances)
        ii = 1.0 / numpy.arange(1,len(idents)+1)
        middle = numpy.sum(offsets * (node_weights*ii)[:,None], 0)/numpy.sum(node_weights*ii)
        axis = numpy.sum((offsets-middle[None,:]) * (travels*node_weights*ii*ii)[:,None], 0)
        axis_length = numpy.sqrt(axis[0]*axis[0]+axis[1]*axis[1])
        if axis_length > 0.0:
            expected_length = numpy.sum(numpy.absolute(travels)*node_weights*ii)
            strength = axis_length / expected_length
            #print axis_length, expected_length
            axis /= axis_length
            offsets -= axis[None,:]*strength * travels[:,None]
            distances *= (1.0-strength)
        
        lengths = numpy.maximum(1e-6, numpy.sqrt(numpy.sum(offsets*offsets, 1)))
        weights = (distances/lengths - 1.0)
        amounts = amount * (1.0-numpy.arange(1,len(idents)+1)/max_n)
    
        this_positions += offsets * (weights*(1.0-numpy.exp(-amounts)))[:,None]    
        self.positions[idents] = this_positions
    
    def iterate(self):
        # Improve positions
        max_ns = float(len(self.links)) / numpy.arange(1,len(self.links)+1)
        random.shuffle(max_ns)
        
        roots = numpy.arange(len(self.links))
        random.shuffle(roots)

        #Working space
        idents = numpy.empty(len(self.links), 'int')
        signs = numpy.empty(len(self.links), 'int')
        distances = numpy.empty(len(self.links), 'float64')
        travels = numpy.empty(len(self.links), 'int')

        for item in self.links:
            random.shuffle(item)

        for i in xrange(len(self.links)):
            self._improve_layout(roots[i], max_ns[i], self.update_amount,  idents,signs,distances,travels)
        
        
        # Lay out nicely
        
        components = [ ]
        used = numpy.zeros(len(self.links), 'bool')

        for i in xrange(len(self.links)):
            if used[i]: continue
            
            component = [ i ]
            used[i] = True
            
            i = 0
            while i < len(component):
                for ident, sign, travel in self.links[component[i]]:
                    if not used[ident]:
                       used[ident] = True
                       component.append(ident)
                i += 1
            
            components.append(numpy.array(component))
            
        components.sort(key=lambda x: len(x), reverse=True)
        
                
        max_y = numpy.sqrt(
                   numpy.maximum.reduce(self.positions[:,0])*
                   numpy.maximum.reduce(self.positions[:,1]) ) * 0.75
        cur_y = 0.0
        left_x = 0.0
        max_x = 0.0
        for component in components:
            cpos = self.positions[component]     
            cpos = cpos - numpy.average(cpos,0)[None,:]
        
            if len(component) > 1:
                U, s, Vh = linalg.svd(cpos, full_matrices=0)
                for i in xrange(len(s)):
                    if numpy.sum(Vh[i,i]) < 0.0: s[i] *= -1        
                cpos = U * s[None,:]
            
            cpos[:,0] -= numpy.minimum.reduce(cpos[:,0])
            cpos[:,1] -= numpy.minimum.reduce(cpos[:,1])
            width = numpy.maximum.reduce(cpos[:,0])
            height = numpy.maximum.reduce(cpos[:,1])
            pad = max(width,1.0) * 0.1
            cpos[:,0] += pad
            cpos[:,1] += pad
            width += pad * 2
            height += pad * 2
            
            if height > max_y: 
                max_y = height
                
            if cur_y+height > max_y:
                left_x = max_x
                cur_y = 0
                        
            cpos[:,0] += left_x
            cpos[:,1] += cur_y
            cur_y += height
            max_x = max(max_x, left_x+width)
            
            self.positions[component] = cpos
        
        self.update_amount *= 1.0 - 1.0/1000

class Graph_viewer:
    def __init__(self, graph, actions=['nothing'], callback=None):
        """
            weights : dictionary mapping name to weight
                      kmers will be colored in rank order of weight
        
        """
        
        self.graph = graph
        self.callback = callback
        
        self.window = gtk.Window()
        self.window.connect('destroy', lambda x: gtk.main_quit())
        self.window.set_default_size(800,600)
        self.window.set_title('Graph viewer')
        
        vbox = gtk.VBox()
        self.window.add(vbox)
        
        self.figure = Figure(figsize=(8,6), dpi=50)
        self.axes = self.figure.add_subplot(111)
        
        
        colors = numpy.empty((len(graph.names), 3))
        sizes = numpy.empty(len(graph.names))
        
        sizes[:] = 2.0
        
        #if weights is None:
        #    #self.axes.plot(graph.positions[:,0], graph.positions[:,1], ',')
        #    
        #    colors[:,:] = [[0.0,0.0,0.0]] 
        #
        #else:
        
        #names = weights.keys()
        #values = weights.values()
        ##names.sort(key=lambda x: weights[x])
        #idents = numpy.array([ graph.name_to_ident[name] for name in names ])
        
        #x = numpy.array(values, dtype='float64')
        x = numpy.array(graph.weights, dtype='float64')
        
        x = numpy.log(x)
        
        x -= numpy.minimum.reduce(x)
        x /= numpy.average(x) * 2.0
        #x /= numpy.sum(x*x)*2.0/numpy.sum(x)
        
        xx = numpy.minimum(x,1.0)
        #x = numpy.arange(len(graph.names)) / float(len(graph.names))
        
        colors[:,0] = 0.5-xx*0.5
        colors[:,1] = 0.75-xx*0.5
        colors[:,2] = 1.0-xx*0.5
        
        sizes[:] = numpy.maximum(x,1.0)**2 #*2.0
            
            #n = 20
            #for i in xrange(n):
            #    start = i*len(names)//n
            #    end = (i+1)*len(names)//n
            #    if start == end: continue
            #    
            #    x = (1.0-float(i)/(n-1)) 
            #    position_block = graph.positions[idents[start:end]]
            #    self.axes.scatter(position_block[:,0],
            #               position_block[:,1],
            #               linewidths=0, 
            #               marker='s',
            #               s=10.0,
            #               c=(0.0,x,x*0.5+0.5),
            #               zorder=i)

        dots = Dots(graph.positions[:,1], graph.positions[:,0], colors, sizes)
        self.axes.add_artist(dots)
        
        #if len(graph.links) < 1000:
        #    for i, (other, other_sign, other_travel) in enumerate(graph.links):
        #        for j in other:
        #            if j > i:
        #                self.axes.plot([graph.positions[i,0],graph.positions[j,0]],
        #                           [graph.positions[i,1],graph.positions[j,1]],
        #                           'k-') 
        
        self.axes.axis('scaled')
        self.axes.set_xlim(0.0, numpy.maximum.reduce(graph.positions[:,0]) * 1.1)
        self.axes.set_ylim(0.0, numpy.maximum.reduce(graph.positions[:,1]) * 1.1)
        
        self.figure.subplots_adjust(top=0.99,bottom=0.05,right=0.99,left=0.05)
        
        #pylab.connect('button_press_event', self._on_click) 
        
        self.annotation_pylab = [ ]
        self.clear_annotation()


        self.canvas = FigureCanvas(self.figure)  # a gtk.DrawingArea
        self.canvas.mpl_connect('button_press_event', self._on_down)
        self.canvas.mpl_connect('button_release_event', self._on_up)
        
        vbox.pack_start(self.canvas)
        
        hbox = gtk.HBox()
        vbox.pack_start(hbox, False, False, 10)
        
        label = gtk.Label('Middle click:')
        hbox.pack_start(label, False, False, 5)
        
        self.radios = { }
        last = None
        for action in actions:
            radio = gtk.RadioButton(group=last, label=action)
            last = radio
            self.radios[action] = radio
            hbox.pack_start(radio, False, False, 5)

        label = gtk.Label('Right click: clear')
        hbox.pack_end(label, False, False, 5)
            
        self.radios[actions[0]].set_active(True)


        toolbar = NavigationToolbar(self.canvas, self.window)
        vbox.pack_start(toolbar, False, False)
        
    def run(self):     
        self.window.show_all()
        gtk.main()        

    def clear_annotation(self):
        self.annotation = { }
      
    def label(self, name, label):
        ident = self.graph.name_to_ident[name]
        self.axes.text(
            self.graph.positions[ident,0],
            self.graph.positions[ident,1],
            label,
            horizontalalignment='center',
            verticalalignment='bottom',
            zorder=100000)

    def arrow(self, names, label):
        positions = [ self.graph.positions[self.graph.name_to_ident[name]] 
                      for name in names if self.graph.has(name) ]
        
        if not positions: return #Error?
        
        max_positions = max(4, (len(positions)+29)//30) #20
        if len(positions) > max_positions:
            positions = [ positions[i*(len(positions)-1)//(max_positions-1)] for i in xrange(max_positions) ]
                
        arrow = Arrow(positions, label, True)
        
        #names = [ name for name in names if self.graph.has(name) ]
        #
        #if len(names) < 2: return #Error?
        #
        #ident1 = self.graph.name_to_ident[names[0]]
        #ident2 = self.graph.name_to_ident[names[-1]]
        #
        #arrow = Arrow(self.graph.positions[ident1],
        #              self.graph.positions[ident2],
        #              label,
        #              True)
        self.axes.add_artist(arrow)                            
        
    def annotate(self, name, mass,r,g,b):
        r *= mass
        g *= mass
        b *= mass
        old_mass, old_r, old_g, old_b = self.annotation.get(name,(0.0,0.0,0.0,0.0))
        self.annotation[name] = (old_mass+mass,old_r+r,old_g+g,old_b+b)

    def refresh_annotation(self):
        while self.annotation_pylab:
            item = self.annotation_pylab.pop(-1)
            item.remove()

        xs = [ ]
        ys = [ ]
        colors = [ ]
        sizes = [ ]
        for name in self.annotation:
            mass,r,g,b = self.annotation[name]
            if not mass: continue
            
            ident = self.graph.name_to_ident[name]
            xs.append(self.graph.positions[ident,0])
            ys.append(self.graph.positions[ident,1])
            colors.append((r/mass,g/mass,b/mass))
            sizes.append(mass)
        
        if xs:
            #thing = self.axes.scatter(
            #    xs,
            #    ys,
            #    s=sizes,
            #    c=colors, 
            #    linewidths=0, 
            #    marker='s',
            #    zorder=10000)
            thing = Dots(numpy.array(ys), numpy.array(xs), numpy.array(colors), numpy.array(sizes), zorder=2)            
            self.axes.add_artist(thing)
            self.annotation_pylab.append(thing)
        
        self.canvas.draw()
    
    def name_from_position(self, x,y):
        xoff = self.graph.positions[:,0] - x
        yoff = self.graph.positions[:,1] - y
        dist2 = xoff*xoff+yoff*yoff
        best = numpy.argmin(dist2)
        
        return self.graph.names[best]
    
    def _on_down(self, event):
        self.down_name = self.name_from_position(event.xdata, event.ydata)
    
    def _on_up(self, event):
        if event.inaxes and event.button == 3:
            self.clear_annotation()
            self.refresh_annotation()
            
        elif event.inaxes and event.button == 2:
            name = self.name_from_position(event.xdata, event.ydata)
            
            if self.callback:
                action = None
                for item in self.radios:
                    if self.radios[item].get_active():
                        action = item

                self.callback(self, action, self.down_name, name)
                self.refresh_annotation()
            
            del self.down_name
            


if __name__ == '__main__':
    n = 50
    graph = Graph(xrange(n))
    for i in xrange(n):
        if random.random() > 0.1:
            graph.link(i, random.randint(n))
    
    for i in xrange(100):
        graph.iterate()

    print graph.positions

    def callback(viewer, action, name):
        print action, name
    
    Graph_viewer(graph, actions=['foo','bar'], callback=callback).run()
