##### Probability function that all elements are zero when n-elements are sampled from m-elements with t non-zero entries  
def Pr(m,t,n):
	temppr = 1
	for i in range(t):
		temppr *= (m-n-i)/(m-i)
	return  temppr ##binomial( m-t ,n)/ binomial(m,n)
##### Cost function for conversion that reduces LPN with m samples to LPN with m-n samples 
def Con(m,t,n):
	return (n^2.8+n^2*(m-n)+ n*(m-n)^2)



##### Probability function that all elements are zero when n-elements are sampled from m-elements under (b, k)-RSD setup
##### where b is the number of blocks, k is the size of each block
def PrforRSD(m,n,b,k):
	l = floor(n/b)
	prob = 1
	for i in range(l):
		prob *= (    1- 1 / (k- i) ) ^b
	for i in range(n- b*l):
		prob *= (1- 1/(k-l))
	return prob


##### parameters
## m : number of samples
## n : dimension of the secret
## t : number of nonzero entries
## thrs : threshold value for each conversion 
## initn : dimension, which will be reduced via the first conversion with thrs =1   
##### Attack complexity with m,t,n,thrs. Here initn is set to reduce the cost in the first main loop.
##### This algorithm outputs an attack complexity and each dimension used in the conversion.




def Quick_RP(m,n,t):
   N = [] 
   Prob = []
   Mu = []
   #for i in range(4):
   N = [ round(n/15*13), round(n/15*1), round(n/15*1)]
   Mu += [N[0]]
   Prob += [1/Pr(m,t,N[0])]
   for i in range(2):
      Mu += [Mu[i]+N[1+i]]
      Prob += [ Prob[i]/Pr(m-Mu[i],t,N[1+i])]
   cost = Con(m,t,N[0])/Pr(m,t,N[0])
   for i in range(1):
      cost += Con(m-Mu[i],t,N[1+i])*Prob[1+i]
   cost += N[2]^2.8*Prob[2]
   mit = RR(log(cost ,2)) 
   return mit 





def Quick_regular_RP(m, n, t):
    b = t
    k = ceil(m/b)
    m = b*k
    Prob = []
    Mu = []
    N = [round(n/15*13), round(n/15*1), round(n/15*1)]
    Mu += [N[0]]
    Prob += [1/Pr_for_regular_LPN(m, N[0], b, k)]
    for i in range(2):
        Mu += [Mu[i] + N[1 + i]]
        Prob += [1/Pr_for_regular_LPN(m, Mu[1 + i], b, k)]
    cost = Con(m, t, N[0]) / Pr_for_regular_LPN(m, N[0], b, k)
    for i in range(2):
        cost += Con(m - Mu[i], t, N[1 + i]) * Prob[1 + i]
    cost += N[2] ^ 2.8 * Prob[2]
    mit = RR(log(cost, 2))
    return mit


