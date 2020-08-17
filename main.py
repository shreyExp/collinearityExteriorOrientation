import numpy as np
import math as mt
import json as js


def doAll():
	cameraVars = {}
	cameraVars["pos"] = getInitialCameraPos()
	cameraVars["angles"] = getInitialCameraAngles()
	pointVars = getPointList()
	cameraVars["rotationMatrix"] = getRotationMatrix(cameraVars["angles"])

def getRotationMatrix(angle):
	omega = angle[0]
	phi = angle[1]
	kappa = angle[2]
	m = np.zeros((3,3))
	m[0,0] = mt.cos(phi) * mt.cos(kappa)
	m[0,1] = mt.sin(omega) * mt.sin(phi) * mt.cos(kappa) + mt.cos(omega)*mt.sin(kappa)
	m[0,2] = -1 * mt.cos(omega) * mt.sin(phi) * mt.cos(kappa) + mt.sin(omega)*mt.sin(kappa)
	m[1,0] = -1 * mt.cos(phi) * mt.sin(kappa)
	m[1,1] = -1 * mt.sin(omega) * mt.sin(phi) * mt.sin(kappa) + mt.cos(omega)*mt.cos(kappa)
	m[1,2] = mt.cos(omega) * mt.sin(phi) * mt.sin(kappa) + mt.sin(omega)*mt.cos(kappa)
	m[2,0] = mt.sin(phi)
	m[2,1] = -1 * mt.sin(omega)*mt.cos(phi)
	m[2,2] = mt.cos(omega)*mt.cos(phi)
	return m

def getRsq(m, rel):
	r = m[0] @ rel
	s = m[1] @ rel
	q = m[2] @ rel
	return(r,s,q)


def getMatrixBee(cameraVars, pointVars):
	noOfPoints = len(pointVars)
	B = np.zeros((noOfPoints*2, 7))
	for i in range(noOfPoints):
		B[2*i] = np.array(blist1(pointVars[i], cameraVars))
		B[2*i + 1] = np.array(blist2(pointVars[i], cameraVars))
	return B

def blist1(pointVar, cameraVars):
	point = pointVar["pos"]
	camera = cameraVars["pos"]
	rel = point - camera
	m = cameraVars["rotationMatrix"]
	f = cameraVars["focalLength"]
	r, s, q = getRsq(m, rel)
	angles = cameraVars["angles"]
	omega = angles[0]
	phi = angles[1]
	kappa = angles[2]
	fq2 = (f/mt.pow(q,2))
	b11 = fq2*(r*(-1*m[2,2]*rel[1] + m[2,1]*rel[2]) - q*(-1*m[0,2]*rel[1] + m[0,1]*rel[2]))

	b12 = fq2*(r*(mt.cos(phi)*rel[0] + mt.sin(omega)*mt.sin(phi)*rel[1] - mt.cos(omega)*mt.sin(phi)*rel[2]))
	b12 = b12  + fq2*(-1*q*(-1*mt.sin(phi)*mt.cos(kappa)*rel[0] + mt.sin(omega)*mt.cos(phi)*mt.cos(kappa)*rel[1]))
	b12 = b12  + fq2*(q*mt.cos(omega)*mt.cos(phi)*mt.cos(kappa)*rel[2])

	b13 = -1*(f/q)*(m[1] @ rel)
	b14 = fq2*(r*m[2,0] - q*m[0,0])

	b15 = fq2*(r*m[2,1] - q*m[0,1])
	b16 = fq2*(r*m[2,2] - q*m[0,2])

	photoPos = pointVar["photoPos"]
	pp = cameraVars["principalPoint"]
	relPhoto = photoPos - pp
	J = relPhoto[0] + f*(r/q)

	l = [b11, b12, b13, -1*b14, -1*b15, -1*b16, J]
	return l

def blist2(pointVar, cameraVars):
	point = pointVar["pos"]
	camera = cameraVars["pos"]
	rel = point - camera
	m = cameraVars["rotationMatrix"]
	f = cameraVars["focalLength"]
	r, s, q = getRsq(m, rel)
	angles = cameraVars["angles"]
	omega = angles[0]
	phi = angles[1]
	kappa = angles[2]
	fq2 = (f/mt.pow(q,2))

	b21 = fq2*(s*(-1*m[2,2]*rel[1] + m[2,1]*rel[2]) - q*(-1*m[1,2]*rel[1] + m[1,1]*rel[2]))

	b22 = fq2*(s*(mt.cos(phi)*rel[0] + mt.sin(omega)*mt.sin(phi)*rel[1] - mt.cos(omega)*mt.sin(phi)*rel[2]))
	b22 = b22 + fq2*(-1*q*(mt.sin(phi)*mt.sin(kappa)*rel[0] - mt.sin(omega)*mt.cos(phi)*mt.sin(kappa)*rel[1]))
	b22 = b22 + fq2*(-1*q*(mt.cos(omega)*mt.cos(phi)*mt.sin(kappa)*rel[2]))

	b23 = (f/q)*(m[0] @ rel)
	b24 = fq2 * (s*m[2,0] - q*m[1,0])
	b25 = fq2 * (s*m[2,1] - q*m[1,1])
	b26 = fq2 * (s*m[2,2] - q*m[1,2])

	photoPos = pointVar["photoPos"]
	pp = cameraVars["principalPoint"]
	relPhoto = photoPos - pp
	K = relPhoto[1] + f*(s/q)
	
	l = [b21, b22, b23, -1*b24, -1*b25, -1*b26, K]
	return l


def Main():
	f = open("points.json", "r")
	s = f.read()
	pointVars = js.loads(s)
	n = len(pointVars)
	for i in range(n):
		pointVars[i]["pos"] = np.array(pointVars[i]["pos"])
		pointVars[i]["photoPos"] = np.array(pointVars[i]["photoPos"])
	cameraVars = getInitialEstimates(pointVars)
	cameraVars["rotationMatrix"] = getRotationMatrix(cameraVars["angles"])
	threshold = 0.01
	thresholdPos = 0.01
	thresholdAngles = 0.001

	print("The initial position of the camera in world coordinates is:\n ")
	print(cameraVars["pos"])
	print("The initial angle of the camera: \n ")
	print(cameraVars["angles"])
	print("\n")
	count = 0
	while(True):
		B = getMatrixBee(cameraVars, pointVars) 
		E = B[:,6]
		B = B[:,0:6]
		delta = np.linalg.lstsq(B, E, rcond=None)[0]
		makeChanges(cameraVars, delta)
		count = count + 1
		print("After " + str(count) + "th iteration")
		print("The delta vector is:")
		print(delta)
		print("The norm of delta is:")
		print(np.linalg.norm(delta))
		print("\n")
		if((np.linalg.norm(delta[:3]) < thresholdAngles) and (np.linalg.norm(delta[3:]) < thresholdPos)):
			break
	print("The final position of the camera in world coordinates is:\n ")
	print(cameraVars["pos"])
	print("The final angle of the camera: \n ")
	print(cameraVars["angles"])


def getInitialEstimates(pointVars):
	n = len(pointVars)
	mat = np.zeros((n*2, 5))
	for i in range(n):
		x = pointVars[i]["photoPos"][0]
		y = pointVars[i]["photoPos"][1]
		X = pointVars[i]["pos"][0]
		Y = pointVars[i]["pos"][1]
		mat[2*i] = np.array([x, -1*y, 1, 0, X])
		mat[2*i + 1] = np.array([y, x, 0, 1, Y])
	col = mat[:,4]
	mat = mat[:,0:4]
	ini = np.linalg.lstsq(mat, col, rcond=None)[0]
	a = ini[0]
	b = ini[1]
	tx = ini[2]
	ty = ini[3]
	omega = 0
	phi = 0
	kappa = mt.atan(a/b)
	fi = open("cameraParameters.json", "r")
	s = fi.read()
	fi.close()
	fdic = js.loads(s)
	cameraVars = {}
	cameraVars["focalLength"] = fdic["focalLength"]
	f = cameraVars["focalLength"]
	cameraVars["principalPoint"] = np.array(fdic["principalPoint"])

	h = 0
	count = 0
	for i in range(n):
		for j in range(i+1, n):
			pho1 = pointVars[i]["photoPos"]
			pho2 = pointVars[j]["photoPos"]
			pd = np.linalg.norm(pho2 -pho1)
			gr1 = pointVars[i]["pos"]
			gr2 = pointVars[j]["pos"]
			gd = np.linalg.norm(gr2 - gr1)
			h = h + f*gd/pd
			count = count + 1
	h = h/count
	cameraVars["pos"] = np.array([tx, ty, h])
	cameraVars["angles"] = np.array([omega, phi, kappa])
	return cameraVars

	

def makeChanges(cameraVars, delta):
	deltaPos = delta[3:]
	deltaAngles = delta[:3]
	cameraVars["pos"] = cameraVars["pos"] + deltaPos
	cameraVars["angles"] = cameraVars["angles"] + deltaAngles
	cameraVars["rotationMatrix"] = getRotationMatrix(cameraVars["angles"])

if __name__ == '__main__':
	Main()
