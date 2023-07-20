
##### Probability function that all elements are zero when n-elements are sampled from m-elements with t non-zero entries  
def Pr(m,t,n):
	return binomial( m-t ,n)/ binomial(m,n)
##### Cost function for conversion that reduces LPN with m samples to LPN with m-n samples 
def Con(m,t,n):
	return (n^2.8+n^2*(m-n)+ n*(m-n)^2)


##### parameters
## m : number of samples
## n : dimension of the secret
## t : number of nonzero entries
## thrs : threshold value for each conversion 
## initn : dimension, which will be reduced via the first conversion with thrs =1   
##### Attack complexity with m,t,n,thrs. Here initn is set to reduce the cost in the first main loop.
##### This algorithm outputs an attack complexity and each dimension used in the conversion.
def Esti(m,t,n,initn,thrs):
	tempset = []
	tempn = initn
	tempsum = initn 
	lvl = 0
	crit = thrs/Pr(m,t,n)
	cost = 0
	while tempsum < n:
		
		while True == 1:
			#tempsum 
			if Con(m-tempsum,t,tempn)/Pr(m,t,tempsum) < crit:
				tempn += 1
				tempsum += 1
			elif Con(m-tempsum,t,tempn)/Pr(m,t,tempsum) > crit:
				tempn -= 1
				tempsum -= 1
				#tempn
				if tempn == 0:
					break
				cost += Con(m-tempsum,t,tempn)/Pr(m,t,tempsum)
				tempset += [tempn]
				lvl += 1
				tempn = 1
				tempsum += 1
		if tempn == 0:
			cost += (n-tempsum)^(2.8)/Pr(m,t,n)
			break
	mit = RR(log(cost ,2)) 
	return mit , tempset


##### Attack complexity with m,t,n. This algorithm finds the optimal thrs with a binary search and outputs a total cost with respect to the thrs. 
def RP(m,n,t):
	initn = 1
	tempsum = 1 
	crit = 1/Pr(m,t,n)
	cost = 0
	while Con(m-tempsum,t,initn)/Pr(m,t,tempsum) < crit: 
		initn += 1
		tempsum += 1
	depth = 10
	stset = [2, round(n )  ,round(n^(1.4))]
	eachcost = [ Esti(m,t,n,initn, stset[0]),Esti(m,t,n,initn, stset[1]), Esti(m,t,n,initn, stset[2])]
	while  depth > 1 :
		cost = min(eachcost)
		if cost == eachcost[0]:
			stset = [stset[0], round( (stset[0]+stset[1])/2 ) ,stset[1]]
			eachcost = [eachcost[0], Esti(m,t,n,initn, stset[1])  ,eachcost[1]]
		elif cost == eachcost[1]:
			stset = [ round( (stset[0]+stset[1])/2 ), stset[1], round( (stset[1]+stset[2])/2 ) ]
			eachcost = [ Esti(m,t,n,initn, stset[0]),  eachcost[1], Esti(m,t,n,initn, stset[2])]
		elif cost == eachcost[2]:
			stset = [stset[1], round( (stset[1]+stset[2])/2 ) ,stset[2]]
			eachcost = [eachcost[1], Esti(m,t,n,initn, stset[1])  ,eachcost[2]]
		depth -= 1
	#print("Cost=", cost)
	return cost






### This algorithm outputs an attack complexity for RSD problem
def RSD(m,n,t):
	return log(Con(m,t,t) + 2^RP(m,n-t,t)[0],2)


















