#aux functions
import numpy as np


def meshgrid2(*arrs):
#this generate array similar to ngrid in matlab for many input vectors
    arrs = tuple(reversed(arrs))  #edit
    lens = map(len, arrs)
    dim = len(arrs)

    sz = 1
    for s in lens:
        sz*=s

    ans = []    
    for i, arr in enumerate(arrs):
        slc = [1]*dim
        slc[i] = lens[i]
        arr2 = np.asarray(arr).reshape(slc) #array Convert the input to an array
        for j, sz in enumerate(lens):
            if j!=i:
                arr2 = arr2.repeat(sz, axis=j) 
        ans.append(arr2)
        
    #print "type ans: ", type(ans) #ans is a list
    i=0
    totallength = np.product(ans[0].shape)
    for variablearray in ans[::-1]:
        ans[i] = np.concatenate(variablearray.reshape(totallength,1))
        i=i+1
    return ans
    
def inplace_change(filename, old_string, new_string):
        s=open(filename).read()
        if old_string in s:
                #print 'Changing "{old_string}" to "{new_string}"'.format(**locals())
                s=s.replace(old_string, new_string)
                f=open(filename, 'w')
                f.write(s)
                f.flush()
                f.close()
        else:
                print 'No occurances of "{old_string}" found.'.format(**locals())
