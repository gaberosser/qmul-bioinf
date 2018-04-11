import collections
import time

# aim: benchmark dictionary access speeds with different methods
# the dictionary values are lists; we need to create a blank list entry when the key doesn't already exist
# we'll just insert toy data each time

# the number of dict keys to add
nkeys = 100000

# the number of data points to add to each dictionary entry
ntoydata = 100

# 1) try ... except
the_dict = dict()

tic = time.time()
for k in xrange(nkeys):
    for t in xrange(ntoydata):
        try:
            the_dict[k].append(t)
        except KeyError:
            the_dict[k] = [t]
toc = time.time()
print "try ... except KeyError: %.3f s" % (toc - tic)

# 2) defaultdict
the_dict = collections.defaultdict(list)
tic = time.time()
for k in xrange(nkeys):
    for t in xrange(ntoydata):
        the_dict[k].append(t)
toc = time.time()
print "collections.defaultdict: %.3f s" % (toc - tic)

# 3) setdefault
the_dict = dict()
tic = time.time()
for k in xrange(nkeys):
    for t in xrange(ntoydata):
        the_dict.setdefault(k, []).append(t)
toc = time.time()
print "setdefault: %.3f s" % (toc - tic)

# 4) if <key> in
the_dict = dict()

tic = time.time()
for k in xrange(nkeys):
    if k not in the_dict:
        the_dict[k] = []
    for t in xrange(ntoydata):
        the_dict[k].append(t)
toc = time.time()
print "if k not in: %.3f s" % (toc - tic)
