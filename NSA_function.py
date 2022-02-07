import pickle
import numpy as np
import xlsxwriter as xl
import matplotlib.pyplot as pyplot
import sys

np.set_printoptions(threshold=sys.maxsize)#np.nan)

def mse(valuesmatrix, groupmatrix, mm, NN, nnv):
    acum=np.zeros(nnv)
    for j in range(mm):
        ind=np.nonzero(groupmatrix==j+1)
        for k in range(nnv):
            emptygroups=0
            aux=valuesmatrix[k,ind[0],ind[1]]
            nk=aux.size
            if(nk>0):
                acum[k]+=aux.var()*nk
            else:
                emptygroups+=1
    return acum/(NN+emptygroups-mm)

def validLocations(groupmatrix,nx,ny):
    validmatrix=np.zeros((ny,nx),dtype=bool)
    ww=0
    for j in range(nx-2):
        for i in range(ny-2):
            valid=np.prod(groupmatrix[i:i+3,j:j+3])
            if valid>0:
                validmatrix[i+1,j+1]=True
                ww+=1
    return validmatrix,ww

def sdvFunc(valuesmatrix, validmatrix, ww, nx, ny):
    acum=0.0
    for j in range(nx-2):
        for i in range(ny-2):
            if validmatrix[i+1,j+1]:
                acum+=valuesmatrix[i:i+3,j:j+3].var()
    if ww==0:
        return 0
    else:
        return acum/float(ww)
        
def calculateOF(valuesmatrix, groupmatrix, mm, NN, nnv, ffdv, rr2, rr2max):
    of=1
    aux=mse(valuesmatrix,groupmatrix,mm,NN,nnv)
    for i in range(nv):
        rr2[i]=1.0-aux[i]/ffdv[i]
        if rr2[i]<0:
            rr2[i]=0
        of*=rr2[i]**rr2max[i]
    return of[0]

def addGroup(valuesmatrix, groupmatrix, validmatrix, mm, NN, nnv, ffdv, rr2, rr2max, nx, ny):
    mof=0.0
    iflag=False
    for j in range(nx-2):
        for i in range(ny-2):
            if(validmatrix[i+1,j+1] and np.prod(groupmatrix[i:i+3,j:j+3]==1)==1):
                aux=groupmatrix[i:i+3,j:j+3].copy()
                groupmatrix[i:i+3,j:j+3]=(mm+1)*np.ones((3,3))
                cof=calculateOF(valuesmatrix, groupmatrix, mm+1, NN, nnv, ffdv, rr2, rr2max)
                groupmatrix[i:i+3,j:j+3]=aux

                if cof>=mof:
                    iflag=True
                    mof=cof
                    mi=i
                    mj=j
    if(iflag):
        return [mof, mi, mj]
    else:
        return [mof]
    
def extendGroup(valuesmatrix, groupmatrix, coordinates, mm, NN, nnv, ffdv, rr2, rr2max, nx, ny):
    mof=0
    iflag=False
    for k in range(mm-1):
        for j in range(len(coordinates)):                #if(groupmatrix[i+1,j+1]==1):
            if(groupmatrix[coordinates[j,0]+1,coordinates[j,1]]==k+2 or groupmatrix[coordinates[j,0],coordinates[j,1]+1]==k+2 or groupmatrix[coordinates[j,0]-1,coordinates[j,1]]==k+2 or groupmatrix[coordinates[j,0],coordinates[j,1]-1]==k+2):
                aux=groupmatrix[coordinates[j,0],coordinates[j,1]]
                groupmatrix[coordinates[j,0],coordinates[j,1]]=k+2
                cof=calculateOF(valuesmatrix, groupmatrix, mm, NN, nnv, ffdv, rr2, rr2max)
                groupmatrix[coordinates[j,0],coordinates[j,1]]=aux

                if cof>=mof:
                    iflag=True
                    mof=cof
                    mi=coordinates[j,0]
                    mj=coordinates[j,1]
                    mk=k+2
    if iflag:
        return [mof, mi, mj, mk]
    else:
        return [mof]

################################################

with open('NSATemp.pickle','rb') as infile:
    zar,ar,z,gdsz,labels,convParams,startingColumn=pickle.load(infile)
nv=z[0]
ngy=z[1]
ngx=z[2]
N=zar[zar!=0].size
m=1
oldof=[]
[var,w]=validLocations(zar, ngx, ngy)
fdv=mse(ar,zar,m,N,nv)
sdv=-1*np.ones(nv)
r2max=-1*np.ones(nv)
for i in range(nv):        
        sdv[i]=sdvFunc(ar[i,:,:],var,w,ngx,ngy)
        r2max[i]=1.0-(sdv[i]/fdv[i])
r2=np.zeros((nv,1))
oldr2=np.zeros((nv,1))
oldof.append(calculateOF(ar, zar, m, N, nv, fdv, r2, r2max))
u=addGroup(ar, zar, var, m, N, nv, fdv, r2, r2max, ngx, ngy)
if(u[0]>oldof[-1]):
    m+=1
    zar[u[1]:u[1]+3,u[2]:u[2]+3]=m*np.ones((3,3))
    oldr2=np.c_[oldr2,r2]
    oldof.append(u[0])
flag=True
while(flag):
    car=np.argwhere(zar[1:-1,1:-1]==1)+1
    u=addGroup(ar, zar, var, m, N, nv, fdv, r2, r2max, ngx, ngy)
    v=extendGroup(ar, zar, car, m, N, nv, fdv, r2, r2max, ngx, ngy)
    if((u[0]-oldof[-1])>8*(v[0]-oldof[-1])):
        if(u[0]>oldof[-1]):
            m+=1
            zar[u[1]:u[1]+3,u[2]:u[2]+3]=m*np.ones((3,3))
            oldr2=np.c_[oldr2,r2]
            oldof.append(u[0])
        else:
            flag=False
    elif(v[0]>oldof[-1]):
        zar[v[1],v[2]]=v[3]
        oldr2=np.c_[oldr2,r2]
        oldof.append(v[0])
    else:
        flag=False
    print(oldof[-1])
################################################

zarMasked=np.ma.masked_where(zar==0, zar)
pyplot.imshow(zarMasked, interpolation='none', origin='lower',extent=[0,ngx*gdsz,0,ngy*gdsz])
pyplot.title('Zones')
pyplot.ylabel('Northing')
pyplot.xlabel('Easting')
pyplot.colorbar()
ax=pyplot.gca()
ax.set_xticks(np.arange(gdsz, gdsz*(ngx+1), gdsz), minor=True)
ax.set_yticks(np.arange(gdsz, gdsz*(ngy+1), gdsz), minor=True)
ax.grid(which='minor', color='k', linestyle='-', linewidth=0.1)
pyplot.savefig('zones.png', dpi=200)
pyplot.close()

pyplot.figure()
pyplot.plot(np.power(oldof,1.0/nv))
pyplot.title('Objective function')
pyplot.ylabel('R^2')
pyplot.xlabel('No. of cells')
pyplot.savefig('OF.png', dpi=200)
pyplot.close()

pyplot.figure()
pyplot.plot(oldr2.T)
pyplot.title('Objective function')
pyplot.legend(labels[startingColumn:])
pyplot.ylabel('R^2')
pyplot.xlabel('No. of cells')
pyplot.savefig('OF2.png', dpi=200)
pyplot.close()

with open('zones.txt', 'w') as outfile:
    zar.tofile(outfile,sep=" ", format="%.5f")

with open('result.pickle', 'wb') as outfile:
    pickle.dump([zar,oldof],outfile)

workbook = xl.Workbook('result.xlsx')
for i in range(nv):
    worksheet=workbook.add_worksheet(labels[i+startingColumn])
    for j in range(ngy):
        worksheet.write_row('A'+str(ngy-j),ar[i,j,:])
worksheet=workbook.add_worksheet('zones')
for j in range(ngy):
    worksheet.write_row('A'+str(ngy-j),zar[j,:])
worksheet=workbook.add_worksheet('stats')
headers=['Minimum','Median','Average','Max','Range','SDV','FDV','R2max']
worksheet.write_row('B1',headers)
ar[ar==0]=np.nan
for i in range(nv):
    worksheet.write_string('A'+str(i+2),labels[i+startingColumn])
    worksheet.write_number('B'+str(i+2),np.nanmin(ar[i,:,:]))
    worksheet.write_number('C'+str(i+2),np.nanmedian(ar[i,:,:]))
    worksheet.write_number('D'+str(i+2),np.nanmean(ar[i,:,:]))
    worksheet.write_number('E'+str(i+2),np.nanmax(ar[i,:,:]))
    worksheet.write_formula('F'+str(i+2),'=E'+str(i+2)+'-B'+str(i+2))
worksheet.write_column('G2',sdv)
worksheet.write_column('H2',fdv)
worksheet.write_column('I2',r2max)
worksheet=workbook.add_worksheet('coordinates')
headers2=['Longitude','Latitude','Zone']
worksheet.write_row('A1',headers2)
for j in range(ngx):
    for i in range(ngy):
        worksheet.write_number('A'+str(j*ngy+i+2),gdsz*j/convParams[0]+convParams[1])
        worksheet.write_number('B'+str(j*ngy+i+2),gdsz*i/convParams[2]+convParams[3])
        worksheet.write_number('C'+str(j*ngy+i+2),zar[i,j])
workbook.close()