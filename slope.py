from shapely.geometry import Point,Polygon,LineString
import matplotlib.pyplot as plt
import numpy as np
import math

class Circle:
    def __init__(self,c,r):
        self.c = np.asarray(c)
        self.r = r
        self.region = Point(c).buffer(r)
        self.outline = self.region.boundary
        self.points = np.array(self.outline.coords)
        
        
class Slope:

    def __init__(self,points):
        points = clock_wise_points(points)
        self.region = Polygon(points)
        self.outline = self.region.boundary
        self.points = np.asarray(self.outline.coords)
        self.NoS = 20
        self.dx = 0.3
        self.dy = 0.3
        self.dr = 0.5
        self.fos = None
        self.minfos = None
        self.slip = None
        self.minslip = None
        self.showed = 0
        

    def set_soil_prop(self,c,fri,uw):
        self.c = c
        self.fri = (fri/180.)*np.pi
        self.uw = uw

    def cut(self,slip):
        N = self.NoS
        self.slip = slip
        block = slip.region.intersection(self.region).boundary
        intersections = np.asarray(slip.outline.intersection(self.outline))
        if len(intersections) != 2:
            return False
        points = self.points
        ctr_x = intersections[:,0]
        ctr_x.sort()
        indices = np.logical_and(points[:,0] > ctr_x[0],points[:,0] < ctr_x[-1])
        inner = points[indices]
        ctr_x = np.hstack((ctr_x,inner[:,0]))
        ctr_x.sort()
        diff = np.diff(ctr_x)
        numbers = (diff/np.sum(diff)*N).astype(int)
        x = np.array([0])
        for i in range(numbers.shape[0]):
            t = np.linspace(ctr_x[i],ctr_x[i+1],numbers[i],endpoint = False)
            x = np.hstack((x,t))
        x = np.hstack((x[1:],ctr_x[-1]))

        n = x.shape[0]
        u_xy = np.zeros((n,2))
        d_xy = np.zeros((n,2))
        u_xy[:,0] = x
        d_xy[:,0] = x
        u_xy[:,1] = self.region.bounds[3] + 0.001
        d_xy[:,1] = self.region.bounds[1] - 0.001
        
        for i in range(1,n-1):
            line = LineString((u_xy[i],d_xy[i]))
            points = np.asarray(line.intersection(block))
            if len(points) != 2:
                return False
            args = np.argsort(points[:,1])
            points = points[args]
            d_xy[i] = points[0]
            u_xy[i] = points[1]
            
        u_xy[0] = intersections[0]
        u_xy[-1] = intersections[-1]
        d_xy[0] = intersections[0]
        d_xy[-1] = intersections[-1]

        area = np.zeros(n-1)
        length = np.zeros(n-1)
        for i in range(n-1):
            p = Polygon((d_xy[i],d_xy[i+1],u_xy[i+1],u_xy[i],d_xy[i]))
            line = LineString((d_xy[i],d_xy[i+1]))
            area[i] = p.area
            length[i] = line.length
            
        diff = np.diff(d_xy,axis = 0)
        alpha = np.arctan2(diff[:,1],diff[:,0])

        self.alpha = alpha
        self.area = area
        self.length = length
        self.u_xy = u_xy
        self.d_xy = d_xy
        return True
    
    def calc_fos(self):
        W = self.uw*self.area
        C = self.c*self.length
        tanfri = math.tan(self.fri)
        fos = (np.sum(C) + np.sum(W*np.cos(self.alpha)*tanfri))\
              /(np.sum(W*np.sin(self.alpha)))
        self.fos = fos
        return fos
    
    def solve(self,search_box):
        xmin,xmax = search_box[0],search_box[2]
        ymin,ymax = search_box[1],search_box[3]
        if ymin < self.region.bounds[3]:
            ymin = self.region.bounds[3] + self.dy
        width = xmax - xmin
        height = ymax - ymin
        nx = int(width/self.dx)
        ny = int(height/self.dy)
        rangex = np.linspace(xmin,xmax,nx)
        rangey = np.linspace(ymin,ymax,ny)
        self.sbx = rangex
        self.sby = rangey
        
        self.minfos =  100
        self.slip_list = []
        for cx in rangex:
            for cy in rangey:
                p = Point(cx,cy)
                rmin = p.distance(self.outline) + self.dr
                bounds = self.region.bounds
                rmax1 = p.distance(Point(bounds[0],bounds[1]))
                rmax2 = p.distance(Point(bounds[2],bounds[3]))
                rmax = min(rmax1,rmax2) - self.dr
                if rmax < rmin:
                    continue
                nr = int((rmax-rmin)/self.dr)
                ranger = np.linspace(rmin,rmax,nr)
                for r in ranger:
                    slip = Circle((cx,cy),r)
                    flag = self.cut(slip)
                    if flag:
                        fos = self.calc_fos()
                        if fos < self.minfos:
                            self.minfos = fos
                            self.minslip = slip
                    self.slip_list.append(slip)
        
    
        
    def show(self,slip = None):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.axis('equal')
        points = self.points
        if slip is None:
            slip = self.minslip
            self.fos = self.minfos
            flag = self.cut(slip)

        flag = self.cut(slip)
        if flag:
            self.calc_fos()
            c = slip.c
            ax.scatter(c[0],c[1],marker = "*",s = 20)
            ax.text(0.7,0.9,"$Fos = %.4f$"%(self.fos),fontsize = 12,transform=ax.transAxes)
            ax.plot(points[:,0],points[:,1],color = "black",lw = 2)
            ax.plot(self.d_xy[:,0],self.d_xy[:,1],color = "red",linestyle = "--")
            ax.plot([c[0],self.d_xy[0][0]],[c[1],self.d_xy[0][1]],c = "pink",lw = 1)
            ax.plot([c[0],self.d_xy[-1][0]],[c[1],self.d_xy[-1][1]],c = "pink",lw = 1)
            for i in range(self.u_xy.shape[0]):
                x = (self.u_xy[i][0],self.d_xy[i][0])
                y = (self.u_xy[i][1],self.d_xy[i][1])
                ax.plot(x,y,linewidth = 0.5,color = "green")
            X, Y = np.meshgrid(self.sbx,self.sby)
            ax.scatter(X,Y,marker = ".",s = 5)
            self.showed += 1
            plt.savefig('%d.jpg'%self.showed)
##            plt.show()
        else:
            print("Slicing  is unsuccessful!")
            
        
        
    
def clock_wise_points(points):
    points = np.asarray(points)
    if points.shape[0] < 3:
        return points
    convex_hull = Polygon(points).convex_hull
    centroid = np.asarray(convex_hull.centroid)
    diff = points - centroid
    rads = np.arctan2(diff[:,0],diff[:,1])
    args = np.argsort(rads)
    return points[args]

if __name__ == "__main__":
    points = [[0,0],[0,5.],[4,5],[6,7],[8,8],[10,10],[15,0],[15,10]]
    c = Circle((6,12),7)
    s = Slope(points)
    s.set_soil_prop(20,10,19)
    s.solve([2,12,5,14])
    s.show()
##    for i in range(len(s.slip_list)):
##        s.show(s.slip_list[i])
    
    
    
