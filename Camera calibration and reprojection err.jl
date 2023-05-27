using Pkg
# pkg.add("MAT")
# pkg.add("ImageEdgeDetection")
# pkg.add("Images")
# pkg.add("Plots")
# pkg.add("ImageFiltering")
# pkg.add("LinearAlgebra")
# pkg.add("Statistics")
# pkg.add("Random")
# pkg.add("Noise")
using Noise 
using Images
using Plots
using ImageFiltering
using Random
using Statistics
using ImageEdgeDetection
using ImageEdgeDetection: Percentile
using LinearAlgebra
using MAT

# We aim to find calibration matrix and reprojection error by perfomig DLT algorithm. Also, Adding normalization to the algorithm to improve the results.

# Part a
data = matread(joinpath(@__DIR__, "cube_points.mat"))
point2d=data["points2d"]
point3d=data["points3d"]
ind=data["connecting_indices"]

p1=scatter(point2d[1,:],point2d[2,:])
p2=scatter(point3d[1,:],point3d[2,:],point3d[3,:])
for i=2:length(ind)
global p2=plot!([point3d[1,ind[i-1]], point3d[1,ind[i]]],[point3d[2,ind[i-1]] ,point3d[2,ind[i]]] ,[point3d[3,ind[i-1]], point3d[3,ind[i]]])
end
display(p1)
display(p2)

# Part b
function calibrate(p2d,p3d)
s=[p3d[1,1] p3d[2,1] p3d[3,1] 1]
A=[s 0 0 0 0 -s.*p2d[1,1];0 0 0 0 s -s.*p2d[2,1]]
for i=2:size(p3d,2)
 s=[p3d[1,i] p3d[2,i] p3d[3,i] 1]
 b=[s 0 0 0 0 -s.*p2d[1,i];0 0 0 0 s -s.*p2d[2,i]]
 A=[A;b]
end
svd_vals=svd(A)
V=svd_vals.V
return V[:,end]
end

N=size(point2d,2)
v=calibrate(point2d,point3d)
M=[v[1:4] v[5:8] v[9:12]]'


# Part c
pnew2=M*[point3d; ones(1,N)]
x=pnew2[1,:]./pnew2[3,:]
y=pnew2[2,:]./pnew2[3,:]

p3=scatter(point2d[1,:],point2d[2,:],mc=:red, ms=5, ma=2)
p3=scatter!(x,y,mc=:blue, ms=5, ma=2)
title!("Obtained 2D points without Normalization")
display(p3)

error=0
for i=1:size(x,2)
 global  error=error+((point2d[1,i]-x[i,1])^2+(point2d[2,i]-y[i,1])^2);
end
error=error/N;
println("Error Value: $error")
# we can see that the error value is small and the points are very close to each other.


# Part d
xhat=(1/N)*sum(point2d[1,:]);
yhat=(1/N)*sum(point2d[2,:]);
dhat=(1/N)*sum(sqrt.((point2d[1,:].-xhat).^2+(point2d[2,:].-yhat).^2));

T=[sqrt(2)/dhat 0 -sqrt(2)*xhat/dhat;0 sqrt(2)/dhat -sqrt(2)*yhat/dhat;0 0 1];
Phat2=T*[(point2d[1,:])';(point2d[2,:])';ones(1,N)]

Xhat=(1/N)*sum(point3d[1,:]);
Yhat=(1/N)*sum(point3d[2,:]);
Zhat=(1/N)*sum(point3d[3,:]);
Dhat=(1/N)*sum(.√((point3d[1,:].-Xhat).^2+(point3d[2,:].-Yhat).^2+(point3d[3,:].-Zhat).^2));

U=[sqrt(3)/Dhat 0 0 -sqrt(3)*Xhat/Dhat;0 sqrt(3)/Dhat 0 -sqrt(3)*Yhat/Dhat;0 0 sqrt(3)/Dhat -sqrt(3)*Zhat/Dhat;0 0 0 1]
Phat3=U*[point3d;ones(1,N)];

v1=calibrate(Phat2[1:2,:],Phat3[1:3,:])
M_normal=[v1[1:4] v1[5:8] v1[9:12]]'
M_unnormal=inv(T)*M_normal*U

p2_unnormal=M_unnormal*[point3d; ones(1,N)]
x1=p2_unnormal[1,:]./p2_unnormal[3,:]
y1=p2_unnormal[2,:]./p2_unnormal[3,:]

p4=scatter(point2d[1,:],point2d[2,:],mc=:red, ms=5, ma=2)
p4=scatter!(x1,y1,mc=:blue, ms=5, ma=2)
title!("Obtained 2D points with Normalization")
display(p4)

error1=0
for i=1:size(x1,2)
 global  error1=error1+((point2d[1,i]-x1[i,1])^2+(point2d[2,i]-y1[i,1])^2);
end
error1=error1/N;
println("Error Value: $error1")
# The error value is slightly smaller than the case we did not use normalization.

