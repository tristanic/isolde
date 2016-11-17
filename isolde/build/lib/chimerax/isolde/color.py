# Define a simple continuous 3-colour scale for colouring objects by some
# property value. Default is red-white-blue. Returns integer arrays defining
# [red, green, blue, opacity].
class ThreeColorScale():
    def __init__(self, minval, maxval, midval = None, mincolor = [255,0,0,255], 
                midcolor = [255,255,255,255], maxcolor = [0,0,255,255]):
        import numpy
        self.minval = minval
        self.maxval = maxval
        if midval is None:
            self.midval = (minval + maxval) / 2
        else:
            self.midval = midval
        
        self.mincolor = numpy.array(mincolor)
        self.midcolor = numpy.array(midcolor)
        self.maxcolor = numpy.array(maxcolor)
    
    # Return the color corresponding to a single value    
    def get_color(self, val):
        if val <= self.minval:
            return self.mincolor
        elif val <= self.midval:
            return ((val-self.minval)/(self.midval-self.minval)*(self.midcolor-self.mincolor)+self.mincolor).astype('uint8')
        elif val <= self.maxval:
            return ((val-self.midval)/(self.maxval-self.midval)*(self.maxcolor-self.midcolor)+self.midcolor).astype('uint8')
        return self.maxcolor
    
    # Return an array of colors corresponding to an array of input values    
    def get_colors(self, vals):
        import numpy    
        c = []
        for v in vals:
            c.append(self.get_color(v))
        return numpy.array(c)

# Some default colour gradients    
color_scales = {
        'RWB': [[255,0,0,255],[255,255,255,255],[0,0,255,255]],
        'BWR': [[0,0,255,255],[255,255,255,255],[255,0,0,255]],
        'RGB': [[255,0,0,255],[0,255,0,255],[0,0,255,255]],
        'BGR': [[0,0,255,255],[0,255,0,255],[255,0,0,255]],
        'RWG': [[255,0,0,255],[255,255,255,255],[0,255,0,255]],
        'GWR': [[0,255,0,255],[255,255,255,255],[255,0,0,255]],      
}

# Returns an object encapsulating one of the standard three-color scales
# defined in color_scales.
def standard_three_color_scale(name, minv, maxv, midv = None):
    minc, midc, maxc = color_scales[name]
    return ThreeColorScale(minv, maxv, midv, 
                minc, midc, maxc)
