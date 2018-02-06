band7=[4,5,6,7]
hasil=[]
for i in band7:
    # a = i < 6
    #hasil.append(a)
    if i > 5:
        hasil.append(i)
    #else:
        #hasil.remove(i)
#print "kosong"
#hasil.append(i)

#hasil.append(i < 6)
# s='print'+ repr(hasil)
#hasil.append(i)


print 'nilai',format(hasil)

import numpy as np
which = lambda lst:list(np.where(lst)[0])
lst = map(lambda band7:band7<5)
print lst