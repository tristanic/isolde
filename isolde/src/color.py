# Copyright 2017 Tristan Croll
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


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
        numpy.seterr(all='raise')
        minc = self.mincolor
        midc = self.midcolor
        maxc = self.maxcolor
        minv = self.minval
        midv = self.midval
        maxv = self.maxval
        invals = numpy.array(vals)
        try:
            above_max = invals >= maxv
            below_min = invals < minv
            between_min_and_mid = numpy.logical_and(invals >=minv, invals < midv)
            between_mid_and_max = numpy.logical_and(invals >=midv, invals < maxv)
        except:
            print(invals)
            raise
        
        
        colors = numpy.empty([len(invals),4], numpy.float32)
        colors[above_max] = maxc
        colors[below_min] = minc
        lm_colors = colors[between_min_and_mid]
        lm_vals = invals[between_min_and_mid]
        mm_colors = colors[between_mid_and_max]
        mm_vals = invals[between_mid_and_max]
        for i in range(4):
            lm_colors[:,i] = (lm_vals-minv)/(midv-minv)*(midc[i]-minc[i])+minc[i]
            mm_colors[:,i] = (mm_vals-midv)/(maxv-midv)*(maxc[i]-midc[i])+midc[i]
        colors[between_min_and_mid] = lm_colors
        colors[between_mid_and_max] = mm_colors
        
        return colors.astype(numpy.uint8)
        

# Some default colour gradients    
COLOR_SCALES = {
        'RWB': [[255,0,0,255],[255,255,255,255],[0,0,255,255]],
        'BWR': [[0,0,255,255],[255,255,255,255],[255,0,0,255]],
        'RGB': [[255,0,0,255],[0,255,0,255],[0,0,255,255]],
        'BGR': [[0,0,255,255],[0,255,0,255],[255,0,0,255]],
        'RWG': [[255,0,0,255],[255,255,255,255],[0,255,0,255]],
        'GWR': [[0,255,0,255],[255,255,255,255],[255,0,0,255]],
        'GYO': [[0,255,0,255],[255,240,50,255],[255,120,50,255]],
        'OYG': [[255,120,50,255],[255,240,50,255],[0,255,0,255]],
        'PiYG': [[255,0,100,255],[255,240,50,255],[0,255,0,255]],
        'GYPi': [[0,255,0,255],[255,240,50,255],[255,0,100,255]],      
}

def standard_three_color_scale(name, minv, maxv, midv = None):
    '''
    Returns an object encapsulating one of the standard three-color scales
    defined in color.COLOR_SCALES.
    Args:
        minv:
            Minimum value - all scores less than this will be coloured the same
        maxv:
            Maximum value - all scores above this will be coloured the same
        midv (optional):
            Mid-point defining where the middle colour sits. If not provided,
            it will be calculated as (minv+maxv)/2
    '''
    minc, midc, maxc = COLOR_SCALES[name]
    return ThreeColorScale(minv, maxv, midv, 
                minc, midc, maxc)
