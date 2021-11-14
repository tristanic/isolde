
class ProblemAggregator:

    def __init__(self, structure):
        self.structure = structure
        self._problem_managers = []
    
    def register_indicator_type(self, problem_getter):
        '''
        Argument "problem_getter" should be a callable taking the structure as its
        sole argument and returning a list of objects representing problem sites. 
        Each object must have a `coord` property returning a single xyz coordinate
        and a `atoms` property returing involved atoms.
        '''
        self._problem_managers.append(problem_getter)
    
    def problem_zones(self, cutoff=3, min_points=6):
        import numpy
        sites = []
        for pm in self._problem_managers:
            sites.extend(pm(self.structure))
        coords = numpy.array([site.coord for site in sites])
        from ..clustering import dbscan
        clusters = dbscan(coords, cutoff, min_points)
        from collections import itemgetter
        clustered_issues = []
        for c in clusters:
            f = itemgetter(*c)
            clustered_issues.append([f(sites), coords[c]])
        return clustered_issues


        




class ProblemAggregatorGUI:

    def __init__(self, isolde):
        pass

