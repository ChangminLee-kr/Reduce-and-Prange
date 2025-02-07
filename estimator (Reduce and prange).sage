

##### Probability function that all elements are zero when n-elements are sampled from m-elements with t non-zero entries  
def Pr(m,t,n):
   return binomial( m-t ,n)/ binomial(m,n)*(1- 1/exp(1))
##### Cost function for conversion that reduces LPN with m samples to LPN with m-n samples 
def Con(m,t,n):
   return (n^2.8+n^2*(m-n)+ n*(m-n)^2)

def Conv(m, n):
   return (n^2.8+n^2*(m-n)+ n*(m-n)^2)



def Pr_for_regular_LPN(m,n,b,k):
	l = floor(n/b)
	prob = 1
	for i in range(l):
		prob *= (    1- 1 / (k- i) ) ^b
	for i in range(n- b*l):
		prob *= (1- 1/(k-l))
	return prob*(1- 1/exp(1))




################### Algebraic RP ###################

def RSD_for_high_degree_F2(m,n,t):
   b = t
   k = round(m/b)
   cmplx = Infinity
   b1 = b - round(b*n/m)
   b0 = round(b*n/m) 
   if b1 == 0:
      for i in range(m-n-1):
         Neq = binomial(m-n-i,2)
         n2 = round(n - (( k-i) +b  + i )) 
         n3 = n - i -n2
         while binomial(m-n - i,2)  > binomial(n3+1,2) - binomial(k -round((i +n2)/b0) ,2)*b0:
               n2 -= 1
               n3 = n - i - n2
         n2 += 1
         n3 = n- i -n2
         p1 = Pr_for_regular_LPN( b0*k  , n2   , b0, k)
         p2 = Pr_for_regular_LPN( m-n  , i , 1, k)
         temp_cmplx = log(Neq^(2.8) + ((i)^2.8 + (i)^2*(m-n - i)+ (i)*(m-n-i)*(n-n2-i))  ,2) - log(p1,2)- log(p2,2)
         print('temp_cmplx=', RR(temp_cmplx) )
         if temp_cmplx < cmplx:
            cmplx = temp_cmplx
   else:
      for i in range(k-1):            #### Guessing n1   k-1 이 맞는데, 편의상 10개로 끊음
         Neq = binomial(k-i,2)*b1
         n2 = round(n - (( k-i)*sqrt(b1) +b  + b1*i ))  ### minimum n2 to guarantee the eqn 1// When F =F2, add b, F is other field remove b
         if n2 < 0:
            break
         n3 = n - i*b1 - n2
         while binomial(k - i,2) * b1 > binomial(n3+1,2) - binomial(k -round((i*b1 +n2)/b0) ,2)*b0:
               n2 -= 1
               n3 = n - i*b1 - n2
         n2 += 1
         n3 = n - i*b1 - n2 
         p1 = Pr_for_regular_LPN( b0*k  , n2   , b0, k)
         p2 = Pr_for_regular_LPN( b1*k  , b1*i , b1, k)
         temp_cmplx = log(Neq^(2.8) + ((b1*i)^2.8 + (b1*i)^2*(m-n - b1*i)+ (b1*i)*(m-n-b1*i)*(n-n2-b1*i))  ,2) - log(p1,2)- log(p2,2)
         print('temp_cmplx=', RR(temp_cmplx) )
         if temp_cmplx < cmplx:
            cmplx = temp_cmplx
         #   print('i=',i)
   return RR(cmplx)


def RSD_for_high_degree(m,n,t):
   b = t
   k = round(m/b)
   cmplx = Infinity
   b1 = b - round(b*n/m)
   b0 = round(b*n/m) 
   if b1 == 0:
      for i in range(m-n-1):
         Neq = binomial(m-n-i,2)
         n2 = round(n - (( k-i) + i )) 
         n3 = n - i -n2
         while binomial(m-n - i,2)  > binomial(n3+1,2) - binomial(k -round((i +n2)/b0) ,2)*b0:
               n2 -= 1
               n3 = n - i - n2
         n2 += 1
         n3 = n- i -n2
         p1 = Pr_for_regular_LPN( b0*k  , n2   , b0, k)
         p2 = Pr_for_regular_LPN( m-n  , i , 1, k)
         temp_cmplx = log(Neq^(2.8) + ((i)^2.8 + (i)^2*(m-n - i)+ (i)*(m-n-i)*(n-n2-i))  ,2) - log(p1,2)- log(p2,2)
         print('temp_cmplx=', RR(temp_cmplx) )
         if temp_cmplx < cmplx:
            cmplx = temp_cmplx
   else:
      for i in range(k-1):            #### Guessing n1   k-1 이 맞는데, 편의상 10개로 끊음
         Neq = binomial(k-i,2)*b1
         n2 = round(n - (( k-i)*sqrt(b1) + b1*i ))  ### minimum n2 to guarantee the eqn 1// When F =F2, add b, F is other field remove b
         if n2 < 0:
            break
         n3 = n - i*b1 - n2
         while binomial(k - i,2) * b1 > binomial(n3+1,2) - binomial(k -round((i*b1 +n2)/b0) ,2)*b0:
               n2 -= 1
               n3 = n - i*b1 - n2
         n2 += 1
         n3 = n - i*b1 - n2 
         p1 = Pr_for_regular_LPN( b0*k  , n2   , b0, k)
         p2 = Pr_for_regular_LPN( b1*k  , b1*i , b1, k)
         temp_cmplx = log(Neq^(2.8) + ((b1*i)^2.8 + (b1*i)^2*(m-n - b1*i)+ (b1*i)*(m-n-b1*i)*(n-n2-b1*i))  ,2) - log(p1,2)- log(p2,2)
         print('temp_cmplx=', RR(temp_cmplx) )
         if temp_cmplx < cmplx:
            cmplx = temp_cmplx
         #   print('i=',i)
   return RR(cmplx)




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
         if Con(m-(tempsum-tempn),t,tempn)/Pr(m,t,tempsum) < crit:
            tempn += 1
            tempsum += 1
            if tempsum == n:
               cost += (Con(m-(tempsum-tempn),t,tempn)/Pr(m,t,tempsum))
               tempset += [tempn]
               tempn = 0
               break   
         elif Con(m-(tempsum-tempn),t,tempn)/Pr(m,t,tempsum) > crit:
            tempn -= 1
            tempsum -= 1
            #tempn
            if tempn == 0:
               break
            cost += (Con(m-(tempsum-tempn),t,tempn)/Pr(m,t,tempsum))
            tempset += [tempn]
            lvl += 1
            tempn = 1
            tempsum += 1
      if tempn == 0:
         cost += ( max( n-tempsum,1)  )^(2.8)/Pr(m,t,n)
         break
   mit = RR(log(cost ,2)) 
   #if len(tempset)>6:
   #   mit = infinity
   return mit , tempset




##### Attack complexity with m,t,n. This algorithm finds the optimal thrs with a binary search and outputs a total cost with respect to the thrs. 
def RP(m,n,t):
   initn = 1
   tempsum = 1 
   crit = 1/Pr(m,t,n)
   cost = 0
   while Con(m,t,initn)/Pr(m,t,tempsum) < crit: 
      initn += 1
      tempsum += 1
   depth = 10
   stset = [2, round(n^1 )  ,round(n^(2.8))]
   eachcost = [ Esti(m,t,n,initn, stset[0]),Esti(m,t,n,initn, stset[1]), Esti(m,t,n,initn, stset[2])]
   while  depth > 1 :
      cost = min(eachcost)
      if cost == eachcost[0]:
         stset = [stset[0], round( sqrt(stset[0]*stset[1]) ) ,stset[1]]
         eachcost = [eachcost[0], Esti(m,t,n,initn, stset[1])  ,eachcost[1]]
      elif cost == eachcost[1]:
         stset = [ round( sqrt(stset[0]*stset[1]) ), stset[1], round( sqrt(stset[1]*stset[2]) ) ]
         eachcost = [ Esti(m,t,n,initn, stset[0]),  eachcost[1], Esti(m,t,n,initn, stset[2])]
      elif cost == eachcost[2]:
         stset = [stset[1], round( sqrt(stset[1]*stset[2]) ) ,stset[2]]
         eachcost = [eachcost[1], Esti(m,t,n,initn, stset[1])  ,eachcost[2]]
      depth -= 1
   #print("Cost=", cost)
   return cost[0]







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





def RP_for_regular_LPN(m,n,t):
	b = t
	k = ceil(m/b)
	m = b*k
	crit = 1/Pr_for_regular_LPN(m,n,b,k)
	cost = 0
	initn = [1, round((1+n)/2), n]
	depth = round(log(n,2)-0.5)
	while depth > 1:
		temp0 = initn[0]
		temp1 = initn[1]
		temp2 = initn[2]
		if Con(m,b,temp1)/Pr_for_regular_LPN(m,temp1,b,k) < crit:
			initn[0] = temp1
			initn[1] = round( (temp2+temp1)/2)
		elif Con(m,b,temp1)/Pr_for_regular_LPN(m,temp1,b,k) > crit:
			initn[2] = temp1
			initn[1] = round( (temp0+temp1)/2)
		depth -= 1
	tempsum = initn[0]
	initn = tempsum
	depth = 15
	stset = [2, round(n )  ,round(n^(2.8))]
	eachcost = [ Esti_for_regular_LPN(m,n,b,k,initn, stset[0]),Esti_for_regular_LPN(m,n,b,k,initn, stset[1]), Esti_for_regular_LPN(m,n,b,k,initn, stset[2])]
	while  depth > 1 :
		cost = min(eachcost)
		if cost == eachcost[0]:
			stset = [stset[0], round( sqrt(stset[0]*stset[1])) ,stset[1]]
			eachcost = [eachcost[0], Esti_for_regular_LPN(m,n,b,k,initn, stset[1])  ,eachcost[1]]
		elif cost == eachcost[1]:
			stset = [ round( sqrt(stset[0]*stset[1]) ), stset[1], round( sqrt(stset[1]*stset[2]) ) ]
			eachcost = [ Esti_for_regular_LPN(m,n,b,k,initn, stset[0]),  eachcost[1], Esti_for_regular_LPN(m,n,b,k,initn, stset[2])]
		elif cost == eachcost[2]:
			stset = [stset[1], round( (stset[1]+stset[2])/2 ) ,stset[2]]
			eachcost = [eachcost[1], Esti_for_regular_LPN(m,n,b,k,initn, stset[1])  ,eachcost[2]]
		depth -= 1
	#print("Cost=", cost)
	return cost[0]



def Esti_for_regular_LPN(m,n,b,k,initn,thrs):
	tempset = []
	tempn = initn
	tempsum = initn 
	lvl = 0
	crit = thrs/Pr_for_regular_LPN(m,n,b,k)
	cost = 0
	while tempsum < n:
		
		while True == 1:
			#tempsum 
			if Con(m-(tempsum-tempn),b,tempn)/Pr_for_regular_LPN(m,tempsum,b,k) < crit:
				tempn += 1
				tempsum += 1
				if tempsum == n:
					cost += (Con(m-(tempsum-tempn),b,tempn)/Pr_for_regular_LPN(m,tempsum,b,k))
					tempset += [tempn]
					tempn = 0
					break	
			elif Con(m-(tempsum-tempn),b,tempn)/Pr_for_regular_LPN(m,tempsum,b,k) > crit:
				tempn -= 1
				tempsum -= 1
				#tempn
				if tempn == 0:
					break
				cost += (Con(m-(tempsum-tempn),b,tempn)/Pr_for_regular_LPN(m,tempsum,b,k))
				tempset += [tempn]
				lvl += 1
				tempn = 1
				tempsum += 1
		if tempn == 0:
			cost += ( max( n-tempsum,1)  )^(2.8)/Pr_for_regular_LPN(m,n,b,k)
			break
	mit = RR(log(cost ,2)) 
	return mit , tempset




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




def RSD_algebraic_RP(m,n,t):
	b = t
	k = round(m/b)
	cmplx = Infinity
	b1 = b - round(b*n/m)
	b0 = round(b*n/m) 
	for i in range(10):				
		Neq = binomial(k-i,2)*b1
		n2 = round(n - (( k-i)*sqrt(b1) +b  + b1*i ))  ### minimum n2 to guarantee the eqn 1// When F =F2, add b, F is other field remove b
		n3 = n - i*b1 - n2
		while binomial(k - i,2) * b1 > binomial(n3+1,2) - binomial(k -round((i*b1 +n2)/b0) ,2)*b0:
				n2 -= 1
				n3 = n - i*b1 - n2
		n2 += 1
		n3 = n - i*b1 - n2 
		p1 = Pr_for_regular_LPN( b0*k  , n2   , b0, k)
		p2 = Pr_for_regular_LPN( b1*k  , b1*i , b1, k)
		temp_cmplx = log(Neq^(2.8) + ((b1*i)^2.8 + (b1*i)^2*(m-n - b1*i)+ (b1*i)*(m-n-b1*i)*(n-n2-b1*i))  ,2) - log(p1,2)- log(p2,2)
		print('temp_cmplx=', RR(temp_cmplx) )
		if temp_cmplx < cmplx:
			cmplx = temp_cmplx
		#	print('i=',i)
	return RR(cmplx)


