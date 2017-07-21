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
        minc = self.mincolor
        midc = self.midcolor
        maxc = self.maxcolor
        minv = self.minval
        midv = self.midval
        maxv = self.maxval
        invals = numpy.array(vals)
        above_max = numpy.argwhere(invals >= maxv).ravel()
        below_min = numpy.argwhere(invals < minv).ravel()
        between_min_and_mid = numpy.argwhere(
            numpy.logical_and(invals >= minv, invals < midv)).ravel()
        between_mid_and_max = numpy.argwhere(
            numpy.logical_and(invals >= midv, invals < maxv)).ravel()
        
        
        outvals = numpy.empty([len(invals),4], numpy.uint8)
        if len(above_max):
            outvals[above_max] = maxc
        if len(below_min):
            outvals[below_min] = minc
        if len(between_min_and_mid):
            lm_vals = invals[between_min_and_mid]
            lm_colors = numpy.empty([len(lm_vals), 4])
            multiplier = (lm_vals - minv)/(midv-minv)
            for i in range(len(lm_vals)):
                lm_colors[i,:] = multiplier[i]*(midc-minc)+minc
            outvals[between_min_and_mid] = lm_colors
        if len(between_mid_and_max):
            mm_vals = invals[between_mid_and_max]
            mm_colors = numpy.empty([len(mm_vals), 4])
            multiplier = (mm_vals-midv)/(maxv-midv)
            for i in range(len(mm_vals)):
                mm_colors[i,:] = multiplier[i]*(maxc-midc)+midc
            outvals[between_mid_and_max] = mm_colors
        
        return outvals

# Some default colour gradients    
color_scales = {
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

# Returns an object encapsulating one of the standard three-color scales
# defined in color_scales.
def standard_three_color_scale(name, minv, maxv, midv = None):
    minc, midc, maxc = color_scales[name]
    return ThreeColorScale(minv, maxv, midv, 
                minc, midc, maxc)
