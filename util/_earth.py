
#from pylab import *
from numpy import array, dot, cross, pi, cos, sin, arccos, arcsin
from numpy.linalg import norm
R_earth = 6371
North = array([0, 0, 1])


def getSC(*args): return getSphericalCoord(*args)
def getSphericalCoord(GC):
    lat = GC.lat / 180. * pi
    long = GC.long / 180. * pi #@ReservedAssignment
    return array([cos(lat) * cos(long), cos(lat) * sin(long), sin(lat)])

def getGC(*args): return getGeographicCoord(*args)
def getGeographicCoord(sCord):
    sCord = sCord / norm(sCord)
    tempLat = arcsin(sCord[2])
    if abs(tempLat) == pi / 2: return 0, 0
    else: return tempLat / pi * 180, arccos(sCord[0] / cos(tempLat)) / pi * 180

class GeographicCoord(object):
    def __init__(self, lat, long): #@ReservedAssignment
        self.lat = lat
        self.long = long

class PlateRotation(object):
    def __init__ (self, rotationCenter, *rotationSpeed):
        """" rotationCenter: GeographicCoord od rotation center of the plate"""
        if len(rotationSpeed) > 0:
            self.rotation = getSphericalCoord(rotationCenter) * rotationSpeed
        else:
            self.rotation = rotationCenter
    def getRotationVector(self):
        return self.rotation
    def getRotationSpeed(self):
        return norm(self.rotation)
    def getRotationCenter(self):
        return getGeographicCoord(self.rotation)

    def getVelocityVectorOfPlateAt(self, POIGeographicCoord):
        """ POI: GeographicCoord of point of interest """
        POI = getSphericalCoord(POIGeographicCoord)
        RC = self.rotation / norm(self.rotation)
        temp = POI - dot(POI, RC) * RC
        distFromRotaionAxis = norm(temp) * R_earth * 1e5
        velocityDirection = cross(RC, POI)
        velocity = velocityDirection / norm(velocityDirection) * norm(self.rotation) * distFromRotaionAxis
        return velocity

def getAzimuthOfVelocity(velVec, POIGeographicCoord):
    POI = getSphericalCoord(POIGeographicCoord)
    NorthAtSurface = cross(POI, cross(North, POI))
    NorthAtSurface = NorthAtSurface / norm(NorthAtSurface)
    return arccos(dot(velVec / norm(velVec), NorthAtSurface)) / pi * 180

def addPlateRotation(A, B):
    return PlateRotation(A.getRotationVector() + B.getRotationVector())




def plates():
    Nazca = PlateRotation(GeographicCoord(56, -94), 7.6e-7 / 180 * pi)
    velVec = Nazca.getVelocityVectorOfPlateAt(GeographicCoord(-28, -71))
    print velVec, norm(velVec)
    print getAzimuthOfVelocity(velVec, GeographicCoord(-28, -71))
    NaPa = PlateRotation(GeographicCoord(55.6, -90.1), 14.2e-7 / 180 * pi)
    PaAnt = PlateRotation(GeographicCoord(-64.3, 96), 9.1e-7 / 180 * pi)
    NaAnt = addPlateRotation(NaPa, PaAnt)
    #print Nazca.getRotationVector()*1e7/pi*180
    #print NaPa.getRotationVector()*1e7/pi*180
    #print PaAnt.getRotationVector()*1e7/pi*180
    #print NaAnt.getRotationVector()*1e7/pi*180
    print NaAnt.getRotationSpeed() * 1e7 / pi * 180, NaAnt.getRotationCenter()

if __name__ == '__main__': # starte Main wenn Modul gestartet wird
    plates()


##Analytische Loesung
#N=1000
#t=linspace(0.0, 2.0*pi, N)
#a0=10
#x0=a0*cos(t)
#z0strich=-a0**2/8.*sin(2*t)
#z0=a0**2/4.*(t-0.5*sin(2*t))
#bx0=-a0*sin(t)
#gamma0=1+0.5*(a0*sin(t))**2
#
#
##Explizites Loesen der PDGL
#t0=t
#N=1000
#t=linspace(0.0, 2.*pi*(1+a0**2/4.), N)
#dt=(max(t)-min(t))/(N-1)
#z=linspace(0., 0., N)
#bz=linspace(0., 0., N)
#ybz=linspace(0., 0., N)
#x=linspace(0., 0., N)
#bx=linspace(0., 0., N)
#ybx=linspace(0., 0., N)
#gamma=linspace(1., 1., N)
#x[0]=a0
#
#for i in arange(N-1):
#    dybx=a0*(-1+bz[i])*cos(t[i]-z[i])*dt 
#    dybz=-a0*bx[i]*cos(t[i]-z[i])*dt     
#    ybx[i+1]=ybx[i]+dybx
#    ybz[i+1]=ybz[i]+dybz
#    gamma[i+1]=(1+ybx[i+1]**2+ybz[i+1]**2)**0.5
#    bx[i+1]=ybx[i+1]/gamma[i+1]
#    bz[i+1]=ybz[i+1]/gamma[i+1]
#    x[i+1]=x[i]+bx[i+1]*dt
#    z[i+1]=z[i]+bz[i+1]*dt
#
#zstrich=z-a0**2/4.*(t-z)
#
#figure(1)
#subplot(411)
#xlabel("z")
#ylabel("x")
#plot(z0strich, x0, label='formel')
#plot(zstrich, x, label='pdgl')
#legend()
#
#subplot(412)
#xlabel("z")
#ylabel("x")
#plot(z0, x0, label='formel')
#plot(z, x, label='pdgl')
#legend()
#
#subplot(413)
#xlabel("t-z")
#ylabel("gammabx")
#plot(t0, bx0, label='formel')
#plot(t-z, gamma*bx, label='pdgl')
#legend()
#
#subplot(414)
#xlabel("t-z")
#ylabel("gamma")
#plot(t0, gamma0, label='formel')
#plot(t-z, gamma, label='pdgl')
#legend()
#
#show()

