import numpy as np
import os
import glob
import matplotlib.pyplot as plt
import copy as cp

###########################
#functions for RK-4 method
def f1(y1,b):
	return y1

def f2(x1,y1,b,sgn,c1):
	if c1 == 0:
		return sgn*(1+y1**2)**(1.5)/(b*(tens0(x1,y1)))#+2*x1)
	else:
		return sgn*(1+y1**2)**(1.5)/(b*(tens1(x1,y1)))
###########################

def tens0(xt,yt): #constant version of surface tension
	return 10

def tens1(xt,yt): #surface tension / surface energy function
	#print('NONCONSTANT!')
	return 10+5*yt

def sec(a,b,atar): #secant method for predicting next 'b' based on current a,b, target_a
	return abs((atar-a[0])*(b[1]-b[0])/(a[1]-a[0])+b[0])

#calculate total interfacial energy of the droplet
def surfenergy(x1,y1,c1,ssl):
	E = 0

	#energy from liquid-vapor component
	for i in range(0,len(x1)-2):
		dx = x1[i+1]-x1[i]
		dy = y1[i+1]-y1[i]
		dl = (dx**2+dy**2)**(0.5)
		if c1 == 0:
			E = E + dl*tens0(x1[i],y1[i])
		else:
			E = E + dl*tens1(x1[i],y1[i])

	#energy from solid-liquid interface at y=0
	#print('liquid surf energy:	' + str(E))
	E = E + ssl*(x1[-1]-x1[0])
	print('total surf energy:	' + str(E))
	return E

#Runge-Kutta 4 algorithm
def RK4(dt,t,z,zp,b,x,sgn,c1):

	k11 = dt*f1(zp,b)
	k12 = dt*f1(zp+1/2*k11,b)
	k13 = dt*f1(zp+1/2*k12,b)
	k14 = dt*f1(zp+k13,b)
	zout = z+1/6*k11+1/3*(k12+k13)+1/6*k14

	k21 = dt*f2(x,zp,b,sgn,c1)
	k22 = dt*f2(x,zp+1/2*k21,b,sgn,c1)
	k23 = dt*f2(x,zp+1/2*k22,b,sgn,c1)
	k24 = dt*f2(x,zp+k23,b,sgn,c1)
	zpout = zp+1/6*k21+1/3*(k22+k23)+1/6*k24

	return(zout,zpout)

def area(x,y): #midpoint integration for area
	A = 0
	for i in range(1,len(x)):
		dx = x[i] -x[i-1]
		A = A + dx*(y[i-1]+y[i])/2
	return A

###########################
#functions for determining which x-y quadrant we're in, for sign of second derivatives
#assumes clockwise rotation through, as a droplet profile generated from -x to +x in +y would follow
def signypp(swapcnt):
	test = swapcnt%4

	if test == 0:
		return 1
	elif test == 1:
		return -1
	elif test == 2:
		return -1
	else:
		return 1

def signxpp(swapcnt):
	test = swapcnt%4

	if test == 0:
		return 1
	elif test == 1:
		return 1
	elif test == 2:
		return -1
	else:
		return -1
###########################

#finds which simulation resulted in smallest energy, SE
def minSE2(SE,l):
	for j in range(1,l):
		if j == 1:
			if (SE[0]) < (SE[1]):
				ind1 = 0
				ind2 = 1
			else:
				ind1 = 1
				ind2 = 0
		else:
			if (SE[j])< (SE[ind1]):
				ind2 = ind1
				ind1 = j
			elif (SE[j]) < (SE[ind2]):
				ind2 = j
	return([ind1,ind2])

def main(xi,xpi,b,n,dl,c1): #RK4 solution method
	#xi = 0 #starting x
	yi = 0 #starting y

	x = np.zeros(n)
	y = np.zeros(n)

	yp = np.zeros(n) #for y(x) frame
	ypp = np.zeros(n)

	xp = np.zeros(n) #for x(y) frame
	xpp = np.zeros(n)

	flag = False
	flagx2 = False
	flagy2 = False
	x[0] = xi
	y[0] = yi

	xp[0] = xpi
	#print(xpi)
	yp[0] = 1/xpi
	xpp[0] = f2(y[0],xp[0],b,-1,c1)

	for i in range(1,n):
		if i == 1:	#start with x(y)
				dy = cp.copy(dl)
				y[i] = y[i-1]+dy
				x[i], xp[i] = RK4(dy,y[i-1],x[i-1],xp[i-1],b,x[i-1],-1,c1)

				sgnxpp = -1
				xpp[i] = f2(y[i],xp[i],b,sgnxpp,c1)

				yp[i] = 1/xp[i]
				sgnypp = np.sign(yp[i])
				ypp[i] = f2(x[i],yp[i],b,sgnypp,c1)
				dx = dl*np.sign(xp[i])

				#determining position in the +/- cycle of innitial condition
				if xp[i]>1:
					flagx = False
					swapcnt = 1
				elif xp[i]<-1:
					flagx = False
					swapcnt = -1
				else:
					flagx = True
					swapcnt = 0
				#print('initial swapcnt')
				#print(swapcnt)
		else:
			if flagx == True:
				#Use x(y) frame
				sgnypp = signypp(swapcnt)
				sgnxpp = signxpp(swapcnt)
				#if i == 2:
				#	print(swapcnt)
				y[i] = y[i-1]+dy
				x[i], xp[i] = RK4(dy,y[i-1],x[i-1],xp[i-1],b,x[i-1],sgnxpp,c1)

				xpp[i] = f2(y[i],xp[i],b,sgnxpp,c1)
				if (flagx2 == False) and (np.sign(xp[i]) != np.sign(xp[i-1])):
					flagx2 = True
					dx = dx*-1
				if abs(xp[i])>1:
					yp[i] = 1/xp[i]
					flagx = False
					flagx2 = False
					swapcnt = swapcnt + 1

			elif flagx == False:
				sgnypp = signypp(swapcnt)
				sgnxpp = signxpp(swapcnt)
				x[i] = x[i-1]+dx
				y[i], yp[i] = RK4(dx,x[i-1],y[i-1],yp[i-1],b,x[i-1],sgnypp,c1)

				ypp[i] = f2(x[i],yp[i],b,sgnypp,c1)

				if (flagy2 == False) and (np.sign(yp[i]) != np.sign(yp[i-1])):
					flagy2 = True
					dy = dy*-1

				if abs(yp[i])>1:
					xp[i] = 1/yp[i]
					flagx = True
					flagy2 = False
					swapcnt = swapcnt + 1

			if (flag == False and y[i] < 0):
#				print(i)
				flag = True
				ifi = i
				break

	if (flag == False):
		ifi = n-1
	return(x,y,ifi)

#routine for converging area of droplet to Atarg, based on modifying initial curvature, b
def Aconv(n,dl,Atarg,c1,xi1,thc1):
	A = [0,0]
	xpi = np.tan((90-thc1)/180*np.pi)	#converting contact angle to dx/dy IC
	print('For contact angle theta_c:	 '+str(thc1))
	print('Goal area:	' + str(Atarg))
	bref = [.1,1] #starting reference radii of curvature
	At = 0
	xt1, yt1, ifi1 = main(xi1,xpi,bref[0],n,dl,c1)
	A[0] = area(xt1[0:ifi1],yt1[0:ifi1])

	xt2, yt2, ifi2 = main(xi1,xpi,bref[1],n,dl,c1)
	A[1] = area(xt2[0:ifi2],yt2[0:ifi2])

	while (abs(At-Atarg)/Atarg) > 1e-5:	#0.001% convergence on contanied area
		bt = sec(A,bref,Atarg)
		xt, yt, ifi = main(xi1,xpi,bt,n,dl,c1)
		At = area(xt[0:ifi],yt[0:ifi])
		print('Current area vs b: ' +str(At) +' vs. ' +str(bt))
		if (abs(A[0] -Atarg)) > (abs(A[1] -Atarg)):
			A[0] = At
			bref[0] = bt
		else:
			A[1] = At
			bref[1] = bt
	return(xt,yt,ifi,bt)

############################################################
## start of initializing & looping
############################################################
#Parameters to change
n 	= 120001 				#maximum number of points
dl 	= 0.0004				#step size in either x or y, depending on frame
Atarg 	= np.pi 				#target constant area
c1 	= 1					# 0 = constant surf tension, otherwise see function tens1()
sigsl 	= 0					#energy value of solid-liquid interface
thcs 	= np.array([40,60,90,120,140])	#left side starting contact angle
xi 	= np.array([0])				#initial conditions for position
nth 	= len(thcs)
#nxi 	= len(xi)

X 	= np.zeros((nth,n))	#x coordinates
Y 	= np.zeros((nth,n))	#y coordinates
IFI 	= np.zeros(nth,dtype=np.dtype(int))	#final used index for coordinates
b 	= np.zeros(nth)	#'b' for secant convergence
SE 	= np.zeros(nth)		#surface energy

for k in range(0,len(xi)):
	thc = np.copy(thcs)
	nth = len(thc)
	print('xi:	' + str(xi[k]))
	for j in range(0,len(thc)): # Solve for each initial angle
		(X[j,:],Y[j,:],IFI[j],b[j]) = Aconv(n,dl,Atarg,c1,xi[k],thc[j]) #converge to desired area from each IC
		SE[j] = surfenergy(X[j][0:IFI[j]],Y[j][0:IFI[j]],c1,sigsl)   #calculate surface energy

	print(SE)
	print(thc)
	print(IFI)
	print(b)

	#index of minimum energy case
	if len(thcs) >= 2:
		minj = minSE2(SE,nth)[0]
	else: minj = 0

	#Plotting converging curves for each initial location + total SE as function of initial angle
	plt.figure(k+1,figsize=(20,10))
	plt.subplot(211)
	for j in range(0,nth):
		if j == minj:
			plt.plot(X[j][0:IFI[j]],Y[j][0:IFI[j]],lw = 3.0, label='$\\theta_c =$'+ "{0:.2f}".format(thc[j]))
		else:
			plt.plot(X[j][0:IFI[j]],Y[j][0:IFI[j]],label='$\\theta_c =$'+ "{0:.2f}".format(thc[j]))

	plt.xlabel('x')
	plt.ylabel('y')
	plt.legend()
	if c1 == 0:
		plt.title('Surface solution family for $\\gamma_{lg} (x) = 10, \\gamma_{ls} = $' + "{0:.2f}".format(sigsl)+', $A = \pi$')
	else:
		plt.title('Surface solution family for $\\gamma_{lg} (x) = 5.0*y+10, \\gamma_{ls} = $' + "{0:.2f}".format(sigsl)+', $A = \pi$')

	plt.grid(True)
	plt.axis('equal')
	plt.plot([-1,8],[0,0],'k')

	plt.subplot(212)
	for j in range(0,nth):
		if j == minj:
			plt.plot(thc[j],SE[j],'o',ms=8.0,label='Surf Energy:'+ "{0:.2f}".format(SE[j]) )
		else:
			plt.plot(thc[j],SE[j],'+',ms=12.0,label='Surf Energy:'+ "{0:.2f}".format(SE[j]) )

	plt.xlabel('$\\theta_c$')
	plt.ylabel('Surface Energy')
	plt.legend()
plt.show()
