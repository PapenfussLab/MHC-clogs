"""
mhc.intervals module

Based largely on an early version of James Taylor's bx.intervals.intersection.py code.
"""


class GenomeIntersector:
    def __init__(self, intersectors=None):
        self.intersectors = {}
        if intersectors:
            self.bulk_add(intersectors)
    
    def add(self, chrom, strand, start, end, name):
        key = (chrom, strand)
        try:
            self.intersectors[key].add(int(start), int(end), value=name)
        except KeyError:
            self.intersectors[key] = Intersector()
            self.intersectors[key].add(int(start), int(end), value=name)
    
    def bulk_add(self, features):
        for f in features:
            self.add(f.chrom, f.strand, f.start, f.end, f.name)
    
    def find(self, chrom, strand, start, end):
        key = (chrom, strand)
        try:
            intersector = self.intersectors[key]
            return self.intersectors[key].find(start, end)
        except KeyError:
            return []
    
    def __getitem__(self, key):
        try:
            intersector = self.intersectors[key]
        except KeyError:
            intersector = Intersector()
        return intersector
    
    def __repr__(self):
        output = []
        for key in self.intersectors:
            for interval in self.intersectors[key].intervals:
                output.append("%s:%i-%i(%s) %s" % (key[0], interval.start, interval.end, key[1], str(interval.value)))
        return "\n".join(output)


class Interval(object):
    """Basic interval, any object with start and end properties will work as well"""
    def __init__(self, start, end, value=None):
        self.start = start
        self.end = end
        self.value = value
    
    def __cmp__(self, other):
        return cmp(self.start, other.start) or cmp(self.end, other.end)
    
    def __repr__(self):
        if self.value:
            return "IntervalWithValue(%d, %d, %s)" % (self.start, self.end, repr(self.value))
        else:
            return "Interval(%d, %d)" % (self.start, self.end)
    
    @staticmethod
    def fromObject(_object, startAttributeName, endAttributeName):
        start = _object.__dict__[startAttributeName]
        end = _object.__dict__[endAttributeName]
        return Interval(start, end, value=_object)


class Intersector(object):
    """
    Data structure for performing window intersect queries on a set of 
    intervals. Algorithm details naively copied from Scott Schwartz's 
    'nxrpts.c'.

    Usage
    =====

    >>> intersecter = Intersector()

    Add intervals, the only requirement is that the interval have integer
    start and end attributes:

    >>> intersecter.add_interval(Interval(0,  10))
    >>> intersecter.add_interval(Interval(3,  7))
    >>> intersecter.add_interval(Interval(3,  40))
    >>> intersecter.add_interval(Interval(10, 50))

    Perform queries:

    >>> intersecter.find(2, 5)
    [Interval(3, 40), Interval(0, 10), Interval(3, 7)]
    >>> intersecter.find(10, 100)
    [Interval(3, 40), Interval(10, 50)]
    >>> intersecter.find(100, 200)
    []
    """

    THRESH = 4 

    # ---- Basic API --------------------------------------------------

    def __init__(self):
        """Initialize"""
        self.intervals = []
        self.dirty = True

    def add_interval(self, interval):
        """Add an interval to the stored set"""
        assert interval.start < interval.end, "Intervals must have length >= 1" 
        self.dirty = True
        self.intervals.append(interval)

    def add(self, start, end, value=None):    
        self.add_interval(Interval(start, end, value))
        
    def find(self, start, end):
        """Return a list of all stored intervals intersecting [start,end)"""
        # Tweak [0,n) to [1,n]
        start += 1
        result = []
        self.find_func(0, len(self.intervals), start, end, result.append)
        return result

    # ---- Preprocessing ----------------------------------------------

    def prepare(self):
        n_intervals = len(self.intervals)
        # starts is filled with the indices of intervals sorted by start position
        self.starts = range(0, n_intervals)
        self.starts.sort(lambda a, b: cmp(self.intervals[a].start+1, self.intervals[b].start+1))
        # ends is filled with the indices of intervals sorted by end position
        self.ends = range(0, n_intervals)
        self.ends.sort(lambda a, b: cmp(self.intervals[a].end, self.intervals[b].end))
        # centers gets filled my partition, initialize it to zeros
        self.centers = [0] * n_intervals
        # Partition recursively
        self.partition(0, n_intervals)
        self.check_partition(0, n_intervals)
        # Ready for queries
        self.dirty = False

    def partition(self, lo, hi):
        if hi - lo < Intersector.THRESH: return
        center = (lo + hi) // 2
        mid = self.intervals[self.ends[center]].end
        # self.ends[lo:center] is correct, separate away self.ends[q:hi)
        q = center + 1
        g = lo
        for i in range(q, hi):
            if self.intervals[self.ends[i]].start+1 > mid:
                self.centers[g] = self.ends[i]
                g += 1
            else:
                self.ends[q] = self.ends[i]
                q += 1
        for i in range(q, hi):
            self.ends[i] = self.centers[lo - q + i]
        # self.starts[q:hi) is correct, separate away self.starts[lo:p)
        p = q
        g = lo
        i = p - 1
        while i >= lo:
            if self.intervals[self.starts[i]].end < mid:
                self.centers[g] = self.starts[i]
                g += 1
            else:
                p -= 1
                self.starts[p] = self.starts[i]
            i -= 1
        for i in range(lo, p):
            self.starts[i] = self.centers[p + lo - 1 - i]
        # Center
        self.centers[center - 1] = p
        self.centers[center] = q
        # Recurse
        self.partition(lo, p)
        self.partition(q, hi) 

    # ---- Find implementation ----------------------------------------

    def find_func(self, lo, hi, start, end, report_func):
        """For each stored interval i that intersects [start, end) call report_func(i)"""
        if self.dirty: self.prepare() 
        if hi - lo < Intersector.THRESH: 
            self.small_find(lo, hi, start, end, report_func)
            return
        p, q, m = self.parts(lo, hi)
        if start > m:
            j = q - 1
            while j >= p and self.intervals[self.ends[j]].end >= start:
                report_func(self.intervals[self.ends[j]])
                j -= 1
            self.find_func(q, hi, start, end, report_func)
        elif end < m:
            j = p
            while j < q and self.intervals[self.starts[j]].start+1 <= end:
                report_func(self.intervals[self.starts[j]])
                j += 1
            self.find_func(lo, p, start, end, report_func)
        else:
            for j in range(p, q):
                report_func(self.intervals[self.ends[j]])
            self.left_find(lo, p, start, end, report_func)
            self.right_find(q, hi, start, end, report_func)
   
    def small_find(self, lo, hi, start, end, report_func):
        for i in range(lo, hi):
            interval = self.intervals[self.starts[i]]
            if start <= interval.end and interval.start+1 <= end:
                report_func(interval)

    def left_find(self, lo, hi, start, end, report_func):
        while hi - lo >= Intersector.THRESH:
            p, q, m = self.parts(lo, hi)
            if start > m:
                j = q - 1
                while j >= p and self.intervals[self.ends[j]].end >= start:
                    report_func(self.intervals[self.ends[j]])
                    j -= 1
                lo = q
            else:
                for j in range(p, hi):
                    report_func(self.intervals[self.ends[j]])
                hi = p
        self.small_find(lo, hi, start, end, report_func)

    def right_find(self, lo, hi, start, end, report_func):
        while hi - lo >= Intersector.THRESH:
            p, q, m = self.parts(lo, hi)
            if end < m:
                j = p
                while j < q and self.intervals[self.starts[j]].start+1 <= end:
                    report_func(self.intervals[self.starts[j]])
                    j += 1
                hi = p
            else:
                for j in range(lo, q):
                    report_func(self.intervals[self.starts[j]])
                lo = q
        self.small_find(lo, hi, start, end, report_func)
            
    def parts(self, lo, hi):
        center = (lo + hi) // 2
        p = self.centers[center - 1]
        q = self.centers[center]
        m = self.intervals[self.ends[center]].end
        return p, q, m

    # ---- Testing ----------------------------------------------------

    def check_partition(self, lo, hi):
        if hi - lo < Intersector.THRESH:
            for i in range(lo, hi - 1):
                assert self.intervals[self.starts[i]].start+1 <= self.intervals[self.starts[i + 1]].start+1
                assert self.intervals[self.ends[i]].end <= self.intervals[self.ends[i + 1]].end
        else:
            p, q, m = self.parts(lo, hi)
            for i in range(lo, p): 
                assert self.intervals[self.ends[i]].end < m
            for i in range(q, hi): 
                assert self.intervals[self.starts[i]].start+1 > m
            for i in range(p, q): 
                assert self.intervals[self.starts[i]].start+1 <= m
                assert self.intervals[self.ends[i]].end >= m
            self.check_partition(lo, p)
            self.check_partition(q, hi)
