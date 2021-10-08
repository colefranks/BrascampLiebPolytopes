import sage.libs.lrcalc.lrcalc as lrcalc
import itertools

#the function BL_vertices([n1,n2,...nt,n]) lists the vertices for 
#generic ni x n B_i's.
#BL_inequalities does the same but for inequalities.

#to test it, run this cell and then the next one.


def add_schurs(list1, list2):
    for x in list2.keys():
        if x in list1.keys():
            list1[x] = list1[x]+list2[x]
        else:
            list1.update({x:list2[x]})

def multiply_schurs(list1, list2, maxrows):
    list3 = {}
    for x in list1.keys():
        for y in list2.keys():
            z = lrcalc.mult(x,y,maxrows)
            add_schurs(list3, z)
    return list3
    
def multiply_many_schurs(poly_list, maxrows):
    p = lrcalc.mult([],[])
    for x in poly_list:
        p = multiply_schurs(p,x,maxrows)
    return p

def dominated(biglist, candidate):
    for listy in biglist:
        dom = true
        for i in range(len(candidate)-1):
            if candidate[i] > listy[i]:
                dom = false
        if candidate[-1] < listy[-1]:
            dom = false
        if dom == true:
            return true
    return false

def numbers_to_partition(bigdim, smalldim, subspace, intersection):
    partition = [bigdim + intersection - smalldim - subspace for i in range(intersection)]
    #partition.extend([0 for i in range(subspace - intersection)])
    return partition
    
def valid(schur,maxcols):
    schurs = [key for key,value in schur.items() if value ==1]
    #print schurs
    if len(schurs)==0:
        return false
    for k in schurs:
        if len(k) == 0:
            return true
    l = [max(partition) for partition in schurs]
    if min(l)<=maxcols:
        return true
    else:
        return false

def is_inequality(bigdim,smalldim_vector, subspace,intersection_vector):
    partitions = [lrcalc.mult([], numbers_to_partition(bigdim, smalldim_vector[i],subspace, intersection_vector[i])) for i in range(len(smalldim_vector))]
    #print partitions
    z = multiply_many_schurs(partitions, subspace)
    #print z
    return valid(z,bigdim - subspace)

def inequalities(dimension_vector):
    dimrange = list(range(1,dimension_vector[-1]))
    allrange = [list(range(i + 1)) for i in dimension_vector[0:-1]]
    for i in allrange:
        i.reverse()
    allrange.append(dimrange)
    potentials = cartesian_product_iterator(allrange)
    inequalities = []
    bigdim = dimension_vector[-1]

    smalldim_vector = dimension_vector[0:-1]

    for potential_inequality in potentials:
        if not dominated(inequalities,potential_inequality):
            #print potential_inequality
            subspace = potential_inequality[-1]
            intersection_vector = potential_inequality[0:-1]
            if is_inequality(bigdim, smalldim_vector, subspace, intersection_vector):
                inequalities.append(potential_inequality)
    return inequalities
    
def positive(m):
    positive = []
    for i in range(m):
        ineq = [0]
        ineq.extend([0 for j in range(i)])
        ineq.append(1)
        ineq.extend([0 for j in range(m - i -1)])
        positive.append(ineq)
    return positive

def BL_vertices(dimensionvector):
    ineqs = inequalities(dimensionvector)
    #print ineqs
    faces = [[ineq[-1]]+[-i for i in ineq[0:-1]] for ineq in ineqs]
    faces.extend(positive(len(dimensionvector)-1))
    p = Polyhedron(ieqs=faces, eqns=[ [dimensionvector[-1]] + [-i for i in dimensionvector[0:-1]] ])
    return p.Vrepresentation()

def P_vertices(dimensionvector):
    ineqs = inequalities(dimensionvector)
    #print(ineqs)
    faces = [[ineq[-1]]+[-i for i in ineq[0:-1]] for ineq in ineqs]
    faces.extend(positive(len(dimensionvector)-1))
    #print(faces)
    full = [dimensionvector[-1]] + [-i for i in dimensionvector[0:-1]]
    faces.append(full)
    #print(faces)
    p = Polyhedron(ieqs=faces)
    return p.Vrepresentation()

#maximizes sum n_i x_i over P, then outputs the face on which it is maximized. 
def P_slice(dimensionvector):
    ineqs = inequalities(dimensionvector)
    #print(ineqs)
    faces = [[ineq[-1]]+[-i for i in ineq[0:-1]] for ineq in ineqs]
    faces.extend(positive(len(dimensionvector)-1))
    #print(faces)
    full = [dimensionvector[-1]] + [-i for i in dimensionvector[0:-1]]
    faces.append(full)
    #print(faces)
    p = Polyhedron(ieqs=faces)
    p.to_linear_program
    lp, x = p.to_linear_program(return_variable=True)
    lp.set_objective(sum(dimensionvector[i]*x[i] for i in range(len(dimensionvector)-1)))
    k = lp.solve()
    print("unweighted maximum on P: ", k)
    newfull = [k] + [-i for i in dimensionvector[0:-1]]
    p = Polyhedron(ieqs=faces, eqns=[newfull])
    #print(lp)
    return p.Vrepresentation()