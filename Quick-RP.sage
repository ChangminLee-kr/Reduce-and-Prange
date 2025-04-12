def Pr_for_regular_LPN(m,n,b,k):
	l = floor(n/b)
	prob = 1
	for i in range(l):
		prob *= (    1- 1 / (k- i) ) ^b
	for i in range(n- b*l):
		prob *= (1- 1/(k-l))
	return prob






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
         n2 = max(0,round(n - (( k-i) +b  + i )) )
         n3 = n - i -n2
         while binomial(m-n - i,2)  > binomial(n3+1,2) - binomial(k -round((i +n2)/b0) ,2)*b0:
               n2 -= 1
               if n2 <0:
                  n2 = 0 
                  break
         n2 += 1
         n3 = n- i -n2
         if n2 < b0*k:
            p1 = Pr_for_regular_LPN( b0*k  , n2   , b0, k)
            p2 = Pr_for_regular_LPN( m-n  , i , 1, k)
            temp_cmplx = log(Neq^(2.8) + m*n + ((i)^2.8 + (i)^2*(m-n - i)+ (i)*(m-n-i)*(n-n2-i))  ,2) - log(p1,2)- log(p2,2) -log( (1-1/e)^2, 2)
            print('temp_cmplx=', RR(temp_cmplx) )
            if temp_cmplx < cmplx:
               cmplx = temp_cmplx
   else:
      for i in ( min(k-1, floor(n/b1) )):            
         Neq = binomial(k-i,2)*b1
         n2 = max(0, round(n - (( k-i)*sqrt(b1) +b  + b1*i )))  ### minimum n2 to guarantee the eqn 1// When F =F2, add b, F is other field remove b
         n3 = n - i*b1 - n2
         while binomial(k - i,2) * b1 > binomial(n3+1,2) - binomial(k -round((i*b1 +n2)/b0) ,2)*b0:
               n2 -= 1
               if n2 < 0:
                n2 = 0
                break
         n2 += 1
         n3 = n - i*b1 - n2 
         if n2 < b0*k:
            p1 = Pr_for_regular_LPN( b0*k  , n2   , b0, k)
            p2 = Pr_for_regular_LPN( b1*k  , b1*i , b1, k)
            temp_cmplx = log(Neq^(2.8) +m*n + ((b1*i)^2.8 + (b1*i)^2*(m-n - b1*i)+ (b1*i)*(m-n-b1*i)*(n-n2-b1*i))  ,2) - log(p1,2)- log(p2,2) -log( (1-1/e)^2, 2)
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
         n2 = max(round(n - (( k-i) + b+ i )) ,0)
         n3 = n - i -n2
         while binomial(m-n - i,2)  > binomial(n3+1,2) - binomial(k -round((i +n2)/b0) ,2)*b0:
               n2 -= 1
               if n2 < 0:
                  n2 =0
                  break
         n2 += 1
         n3 = n- i -n2
         if n2 < b0*k:
            p1 = Pr_for_regular_LPN( b0*k  , n2   , b0, k)
            p2 = Pr_for_regular_LPN( m-n  , i , 1, k)
            temp_cmplx = log(Neq^(2.8) + m*n + ((i)^2.8 + (i)^2*(m-n - i)+ (i)*(m-n-i)*(n-n2-i))  ,2) - log(p1,2)- log(p2,2) -log( (1-1/e)^2, 2)
            print('temp_cmplx=', RR(temp_cmplx) )
            if temp_cmplx < cmplx:
               cmplx = temp_cmplx
   else:
      for i in range( min(k-1, floor(n/b1) ))  :          
         Neq = binomial(k-i,2)*b1
         n2 = max(0, round(n- (k-i)*sqrt(b1) + b1*i + b  ))  ### minimum n2 to guarantee the eqn 1// When F =F2, add b, F is other field remove b
         n3 = max(0, n-i*b1 -n2)
         while binomial(k - i,2) * b1 > binomial(n3+1,2) - binomial(k -round((i*b1 +n2)/b0) ,2)*b0:
               n2 -= 1
               n3 = n - i*b1 - n2
               if n2 < 0:
                  n2 = 0
                  break
         n2 += 1
         n3 = n - i*b1 - n2 
         if n2 < b0*k:
            p1 = Pr_for_regular_LPN( b0*k  , n2   , b0, k)
            p2 = Pr_for_regular_LPN( b1*k  , b1*i , b1, k)
            temp_cmplx = log(Neq^(2.8) +m*n + ((b1*i)^2.8 + (b1*i)^2*(m-n - b1*i)+ (b1*i)*(m-n-b1*i)*(n-n2-b1*i))  ,2) - log(p1,2)- log(p2,2) -log( (1-1/e)^2, 2)
            print('temp_cmplx=', RR(temp_cmplx) )
            if temp_cmplx < cmplx:
               cmplx = temp_cmplx
            #   print('i=',i)
   return RR(cmplx)
