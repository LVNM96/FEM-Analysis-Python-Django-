from django.shortcuts import render
from django.views import View
from .forms import XYZForm, BCForm, FORCEForm, ELDATAForm, UploadForm
import numpy as np
import matplotlib.pyplot as plt
from django.http import HttpResponse
from django.core.files.storage import default_storage        
 

class Index(View):
    def get(self,request):
        formXYZ=XYZForm()
        formBC=BCForm()
        formFORCE=FORCEForm()
        formELDATA=ELDATAForm()
        formUpload = UploadForm()
        return render(request, 'prvi/LS_index.html', {'formXYZ': formXYZ, 'formBC':formBC, 'formFORCE':formFORCE, 'formELDATA':formELDATA, 'formUpload': formUpload,})

    def post(self,request):
        formXYZ=XYZForm(request.POST)
        formBC=BCForm(request.POST)
        formFORCE=FORCEForm(request.POST)
        formELDATA=ELDATAForm(request.POST)
        formUpload = UploadForm(request.POST, request.FILES)
                
        if formUpload.is_valid():
            
            upload = formUpload.cleaned_data['file']
            
            sadrzaj = upload.read().decode('utf-8')

            lines = sadrzaj.split('\n')
            def el_te2d(ie, ELDATA, x1, x2, x3, y1, y2, y3):
                
                YE = ELDATA[0]
                nu = ELDATA[1]
                Debljina = ELDATA[2]            
                
                Det = (x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1)

                if ( Det > 0.000001):
                    pass
                else:
                    ERR_ELEM = 1
                    print(ERR_ELEM)
                
                alpha1 = (x2*y3-x3*y2)/Det
                beta1  = -(y3-y2)/Det
                gama1  = (x3-x2)/Det
                alpha2 = (x3*y1-x1*y3)/Det
                beta2  = -(y1-y3)/Det
                gama2  = (x1-x3)/Det
                alpha3 = (x1*y2-x2*y1)/Det
                beta3  = -(y2-y1)/Det
                gama3  = (x2-x1)/Det
                        
                Be = np.array([[beta1, 0, beta2, 0, beta3, 0],
                    [0, gama1, 0, gama2, 0, gama3],
                    [gama1, beta1, gama2, beta2, gama3, beta3]])
                DeInv = 1/YE * np.array([[1, -nu, 0] , 
                    [-nu, 1, 0], 
                    [0, 0, 2*(1+nu)]])
                De=np.linalg.inv(DeInv)
                Ve = Det * Debljina / 2
                k  = np.dot(np.dot(np.transpose(Be),De),np.dot(Be,Ve)) 
                return k, De, Be
            global DISP    
            
            global KG

            global Stress
            
            XYZ_vrijednosti = [list(map(int, line.split())) for line in lines[:4]]
            XYZ = np.array(XYZ_vrijednosti)
            
            BC_vrijednosti = [list(map(int, line.split())) for line in lines[4:8]]
            BC = np.array(BC_vrijednosti)      


            for broj in BC.flatten():
                if broj != 0 and broj != 1:
                    greska1 = "Greška: Pogrešno unesene vrijednosti matrice BC (uvjeti oslanjanja)"
                    return render(request, 'prvi/error.html', {'greska1':greska1})

            eldata_values = list(map(float, lines[12].split()))
            ELDATA = np.array(eldata_values)   
            
            if ELDATA[0]<0:
                greska2="Greška: Unesene negativne vrijednosti Youngovog modula elastičnosti"
                return render(request, 'prvi/error.html', {'greska2':greska2})  
            elif ELDATA[1]<0:
                greska3="Greška: Unesene negativne vrijednosti Poissonovog koeficijenta"
                return render(request, 'prvi/error.html', {'greska3':greska3})          
            elif ELDATA[2]<0:
                greska4="Greška: Unesene negativne vrijednosti debljine"
                return render(request, 'prvi/error.html', {'greska4':greska4})  
             
            NODSPE = np.array([[0, 1, 2],     #numeriranje elemenata
            [0, 2, 3]])            

            NODS = np.size(XYZ, 0)   #broj cvorova jednak broju redaka XYZ
            NDOF=np.size(BC,1)    #broj stupnjeva slobode jednak broju stupaca BC
            NDIM = np.size(XYZ, 1)  #dimenzija prostora jednaka broju stupaca XYZ

            FORCE = np.zeros((NODS,NDOF))   #za sile zadajemo matricu sa nulama za pocetak
                       
            FORCE_vrijednosti = [list(map(int, line.split())) for line in lines[8:12]]
            FORCE = np.array(FORCE_vrijednosti)          
            
            L=np.amax(XYZ[:,0])   # duljina grede i sirina grede   
            H=np.amax(XYZ[:,1])
           
            DISP = np.zeros((NODS,NDOF)) # pomaci

            NDOFEN = 2  # broj stupnjeva slobode po cvoru

            NEL = np.size(NODSPE,0)        # broj elemenata
            NODE = np.size(NODSPE,1)       # broj čvorova po elementu 

            NEQ=0
            for i in range (0,NODS): 
                for j in range (0,NDOF):
                    if (BC[i][j]==0 ): 
                        NEQ+=1
                        BC[i][j]=NEQ
                    elif (BC[i][j]==1):
                        BC[i,j]=-1
                    else:
                        ERR=1          

            OSL = 0
            for i in range(NODS):
                for j in range(NDOF):
                    if BC[i][j] == -1:
                        OSL += 1
                        BC[i][j] = -(NEQ + OSL)
                    elif BC[i][j] == 0:
                        ERR = 2
                     
            GNEQ = NEQ + OSL 

            KG = np.zeros((GNEQ, GNEQ)) # globalna matrica krutosti 

            NDOFE = NODE*NDOFEN

            for ie in range(NEL): # for petlja za izracun globalne matrice krutosti
    
                x1 = XYZ[NODSPE[ie,0], 0]
                x2 = XYZ[NODSPE[ie,1], 0]
                x3 = XYZ[NODSPE[ie,2], 0]
                y1 = XYZ[NODSPE[ie,0], 1]
                y2 = XYZ[NODSPE[ie,1], 1]
                y3 = XYZ[NODSPE[ie,2], 1]
                k, De, Be = el_te2d(ie, ELDATA, x1,x2,x3,y1,y2,y3)
                
                n=0
                LM=[]
                
                for i in range(NODE):  
                    for j in range(NDOFEN):
                        n=n+1
                        LM.append(BC[NODSPE[ie,i],j])
                        
                    
                for i in range(NDOFE): 
                    IG = np.abs(LM[i]) 
                    for j in range(NDOFE):
                        JG = np.abs(LM[j])
                        if (IG > 0 and JG > 0):
                            
                            
                            KG[IG-1][JG-1] = KG[IG-1][JG-1] + k[i][j]
                        else:
                            ERR = 3
            KGJJ = KG[NEQ:NEQ+OSL, NEQ:NEQ+OSL] 
            KGII = KG[0:NEQ,0:NEQ]
            KGIJ = KG[0:NEQ, NEQ:NEQ+OSL]
            KGJI = KGIJ.transpose()

            # matrica sila
            F=np.zeros((NEQ+OSL,1))
            for i in range(NODS):
                for j in range(NDOF):
                    k = np.abs(BC[i][j]) 
                    F[k-1][0]= FORCE[i][j]

            FI = F[:NEQ,:]

            # matrica pomaka
            D = np.zeros((NEQ+OSL, 1))
            for i in range(NODS):
                for j in range(NDOF):
                    k = np.abs(BC[i][j])   
                    
                    
                    D[k-1][0]= DISP[i][j] 

            DJ = D[NEQ:NEQ+OSL, :]

            DI =np.linalg.inv(KGII) @ (FI - KGIJ @ DJ)

            FJ = (KGJI@DI) + (KGJJ@DJ)

            D = np.vstack((DI, DJ))

            F = np.vstack((FI, FJ))

            for i in range(NODS):
                for j in range(NDOF):
                    k = np.abs(BC[i][j]) 
                    
                    FORCE[i][j] = F[k-1]

            for i in range(NODS):
                for j in range(NDOF):
                    k = abs (BC[i,j])    
                    DISP[i][j] = D[k-1]

            Stress = np.zeros((3, NEL))
            sigma_Mises=np.zeros((2))
            for ie in range(NEL): # izracun deformacija i naprezanja
                n=0
                LM=[]
                for i in range (NODE): 
                    for j in range (NDOFEN):
                        
                        n = n+1
                        
                        LM.append(BC[NODSPE[ie,i],j])
                        
                        ng = np.abs(LM)
                        
                        dispE=[n,0]
                        dispE=D[ng-1]                   
              
            
                x1 = XYZ[NODSPE[ie,0], 0]
                x2 = XYZ[NODSPE[ie,1], 0]
                x3 = XYZ[NODSPE[ie,2], 0]
                y1 = XYZ[NODSPE[ie,0], 1]
                y2 = XYZ[NODSPE[ie,1], 1]
                y3 = XYZ[NODSPE[ie,2], 1]
                [ k, De, Be ] = el_te2d( ie, ELDATA, x1,x2,x3,y1,y2,y3 )
                global sigmaE
                epsE = Be @ dispE
                
                sigmaE = De @ epsE  
                
                x_koord=XYZ[:,0]
                y_koord=XYZ[:,1]
                x_koord=np.append(x_koord,x_koord[0])
                y_koord=np.append(y_koord,y_koord[0]) 
                DISP2=10000*np.array(DISP)
                XYZ_2=np.add(XYZ,DISP2)
                x_koord2=XYZ_2[:,0]
                y_koord2=XYZ_2[:,1]
                x_koord2=np.append(x_koord2,x_koord2[0])
                y_koord2=np.append(y_koord2,y_koord2[0])
                plt.plot(x_koord, y_koord, 'b', linewidth=0.5)
                plt.plot(x_koord2, y_koord2, 'r', linewidth=0.5)
                plt.xlabel('X-os')
                plt.ylabel('Y-os')
                plt.title('Ploča nakon deformacije')
                plt.legend(['početno', 'krajnje(pomaci*10000)'], bbox_to_anchor=(1.15, 1.15))
                plt.grid(True)
                plt.savefig('prvi/static/prvi/xyz_plot.png')
                
                Stress[:, ie] = sigmaE.flatten()         
               
                
                sigma_Mises[ie]= np.sqrt(sigmaE[0]**2 + sigmaE[1]**2 - sigmaE[0]*sigmaE[1] + 3*sigmaE[2]**2)
            varijable={
                'XYZ':XYZ,
                'KG':KG,
                'FORCE':FORCE,
                'DISP':DISP,
                'BC':BC,
                'ELDATA':ELDATA,
                'sigmaE':sigmaE,
                'sigma_Mises':sigma_Mises,
                'xyz_plot_path': '/prvi/static/prvi/xyz_plot.png', 
                'Stress':Stress,
                
            }
            return render(request, 'prvi/LS_rezultati.html', varijable)
 
        elif formXYZ.is_valid() and formBC.is_valid() and formFORCE.is_valid() and formELDATA.is_valid():
            
            
            def el_te2d(ie, ELDATA, x1, x2, x3, y1, y2, y3):
                
                YE = ELDATA[0]
                nu = ELDATA[1]
                Debljina = ELDATA[2]              
             
                Det = (x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1)

                if ( Det > 0.000001):
                    pass
                else:
                    ERR_ELEM = 1
                    print(ERR_ELEM)
                

                alpha1 = (x2*y3-x3*y2)/Det
                beta1  = -(y3-y2)/Det
                gama1  = (x3-x2)/Det
                alpha2 = (x3*y1-x1*y3)/Det
                beta2  = -(y1-y3)/Det
                gama2  = (x1-x3)/Det
                alpha3 = (x1*y2-x2*y1)/Det
                beta3  = -(y2-y1)/Det
                gama3  = (x2-x1)/Det
                        
                Be = np.array([[beta1, 0, beta2, 0, beta3, 0],
                    [0, gama1, 0, gama2, 0, gama3],
                    [gama1, beta1, gama2, beta2, gama3, beta3]])
                DeInv = 1/YE * np.array([[1, -nu, 0] , 
                    [-nu, 1, 0], 
                    [0, 0, 2*(1+nu)]])
                De=np.linalg.inv(DeInv)
                Ve = Det * Debljina / 2
                k  = np.dot(np.dot(np.transpose(Be),De),np.dot(Be,Ve)) 
                return k, De, Be

            
            XYZ = np.array([[formXYZ.cleaned_data['x1'], formXYZ.cleaned_data['y1']],
                            [formXYZ.cleaned_data['x2'], formXYZ.cleaned_data['y2']],
                            [formXYZ.cleaned_data['x3'], formXYZ.cleaned_data['y3']],
                            [formXYZ.cleaned_data['x4'], formXYZ.cleaned_data['y4']]])         
                       
            BC = np.array([[formBC.cleaned_data['BCx1'], formBC.cleaned_data['BCy1']],
                            [formBC.cleaned_data['BCx2'], formBC.cleaned_data['BCy2']],
                            [formBC.cleaned_data['BCx3'], formBC.cleaned_data['BCy3']],
                            [formBC.cleaned_data['BCx4'], formBC.cleaned_data['BCy4']]])
            
            for broj in BC.flatten():
                if broj != 0 and broj != 1:
                    greska1 = "Greska: Pogresno unesene vrijednosti matrice BC"
                    return render(request, 'prvi/error.html', {'greska1':greska1})

            
            
            NODSPE = np.array([[0, 1, 2],     #numeriranje elemenata
            [0, 2, 3]])
         
            NODS = np.size(XYZ, 0)   #broj cvorova jednak broju redaka XYZ
            NDOF=np.size(BC,1)    #broj stupnjeva slobode jednak broju stupaca BC
            NDIM = np.size(XYZ, 1)  #dimenzija prostora jednaka broju stupaca XYZ

            FORCE = np.zeros((NODS,NDOF))   #za sile zadajemo matricu sa nulama za pocetak
            FORCE = np.array([[formFORCE.cleaned_data['Fx1'], formFORCE.cleaned_data['Fy1']],
                            [formFORCE.cleaned_data['Fx2'], formFORCE.cleaned_data['Fy2']],
                            [formFORCE.cleaned_data['Fx3'], formFORCE.cleaned_data['Fy3']],
                            [formFORCE.cleaned_data['Fx4'], formFORCE.cleaned_data['Fy4']]])
            L=np.amax(XYZ[:,0])   # duljina grede i sirina grede   
            H=np.amax(XYZ[:,1])
            
           


            DISP = np.zeros((NODS,NDOF)) # pomaci

            NDOFEN = 2  # broj stupnjeva slobode po cvoru

            NEL = np.size(NODSPE,0)        # broj elemenata
            NODE = np.size(NODSPE,1)       # broj čvorova po elementu 
   
            NEQ=0
            for i in range (0,NODS): 
                for j in range (0,NDOF):
                    if (BC[i][j]==0 ): 
                        NEQ+=1
                        BC[i][j]=NEQ
                    elif (BC[i][j]==1):
                        BC[i,j]=-1
                    else:
                        ERR=1
        
            OSL = 0
            for i in range(NODS):
                for j in range(NDOF):
                    if BC[i][j] == -1:
                        OSL += 1
                        BC[i][j] = -(NEQ + OSL)
                    elif BC[i][j] == 0:
                        ERR = 2
                    
            ELDATA = []
            for ii in range(NEL):
                ELDATA = [formELDATA.cleaned_data['YE'], formELDATA.cleaned_data['nu'], formELDATA.cleaned_data['Debljina']] 

            if ELDATA[0]<0:
                greska2="Unesene negativne vrijednosti Youngovog modula elastičnosti"
                return render(request, 'prvi/error.html', {'greska2':greska2})  
            elif ELDATA[1]<0:
                greska3="Unesene negativne vrijednosti Poissonovog koeficijenta"
                return render(request, 'prvi/error.html', {'greska3':greska3})          
            elif ELDATA[2]<0:
                greska4="Unesene negativne vrijednosti debljine"
                return render(request, 'prvi/error.html', {'greska4':greska4})  
            GNEQ = NEQ + OSL 

            KG = np.zeros((GNEQ, GNEQ)) # globalna matrica krutosti 

            NDOFE = NODE*NDOFEN

            for ie in range(NEL): # for petlja za izracun globalne matrice krutosti
    
                x1 = XYZ[NODSPE[ie,0], 0]
                x2 = XYZ[NODSPE[ie,1], 0]
                x3 = XYZ[NODSPE[ie,2], 0]
                y1 = XYZ[NODSPE[ie,0], 1]
                y2 = XYZ[NODSPE[ie,1], 1]
                y3 = XYZ[NODSPE[ie,2], 1]
                k, De, Be = el_te2d(ie, ELDATA, x1,x2,x3,y1,y2,y3)
                

                n=0
                LM=[]
                
                for i in range(NODE):  
                    for j in range(NDOFEN):
                        n=n+1
                        LM.append(BC[NODSPE[ie,i],j])
                        
                    
                for i in range(NDOFE): 
                    IG = np.abs(LM[i]) 
                    for j in range(NDOFE):
                        JG = np.abs(LM[j])
                        if (IG > 0 and JG > 0):
                            
                            
                            KG[IG-1][JG-1] = KG[IG-1][JG-1] + k[i][j]
                        else:
                            ERR = 3
            KGJJ = KG[NEQ:NEQ+OSL, NEQ:NEQ+OSL] 
            KGII = KG[0:NEQ,0:NEQ]
            KGIJ = KG[0:NEQ, NEQ:NEQ+OSL]
            KGJI = KGIJ.transpose()

            # matrica sila
            F=np.zeros((NEQ+OSL,1))
            for i in range(NODS):
                for j in range(NDOF):
                    k = np.abs(BC[i][j]) 
                    F[k-1][0]= FORCE[i][j]

            FI = F[:NEQ,:]

            # matrica pomaka
            D = np.zeros((NEQ+OSL, 1))
            for i in range(NODS):
                for j in range(NDOF):
                    k = np.abs(BC[i][j])   
                    
                    
                    D[k-1][0]= DISP[i][j] 

            DJ = D[NEQ:NEQ+OSL, :]

            DI =np.linalg.inv(KGII) @ (FI - KGIJ @ DJ)

            FJ = (KGJI@DI) + (KGJJ@DJ)

            D = np.vstack((DI, DJ))

            F = np.vstack((FI, FJ))

            for i in range(NODS):
                for j in range(NDOF):
                    k = np.abs(BC[i][j]) 
                    
                    FORCE[i][j] = F[k-1]

            for i in range(NODS):
                for j in range(NDOF):
                    k = abs (BC[i,j])    
                    DISP[i][j] = D[k-1]

            Stress = np.zeros((3, NEL))
            sigma_Mises=np.zeros((2))
            for ie in range(NEL): # izracun deformacija i naprezanja
                n=0
                LM=[]
                for i in range (NODE): 
                    for j in range (NDOFEN):
                        
                        n = n+1
                        
                        LM.append(BC[NODSPE[ie,i],j])
                        
                        ng = np.abs(LM)
                        
                        dispE=[n,0]
                        dispE=D[ng-1]
                        
                                
                x1 = XYZ[NODSPE[ie,0], 0]
                x2 = XYZ[NODSPE[ie,1], 0]
                x3 = XYZ[NODSPE[ie,2], 0]
                y1 = XYZ[NODSPE[ie,0], 1]
                y2 = XYZ[NODSPE[ie,1], 1]
                y3 = XYZ[NODSPE[ie,2], 1]
                [ k, De, Be ] = el_te2d( ie, ELDATA, x1,x2,x3,y1,y2,y3 )
                
                epsE = Be @ dispE
                
                sigmaE = De @ epsE  
                
                x_koord=XYZ[:,0]
                y_koord=XYZ[:,1]
                x_koord=np.append(x_koord,x_koord[0])
                y_koord=np.append(y_koord,y_koord[0]) 
                DISP2=10000*np.array(DISP)
                XYZ_2=np.add(XYZ,DISP2)
                x_koord2=XYZ_2[:,0]
                y_koord2=XYZ_2[:,1]
                x_koord2=np.append(x_koord2,x_koord2[0])
                y_koord2=np.append(y_koord2,y_koord2[0])
                plt.plot(x_koord, y_koord, 'b', linewidth=0.5)
                plt.plot(x_koord2, y_koord2, 'r', linewidth=0.5)
                plt.xlabel('X-os')
                plt.ylabel('Y-os')
                plt.title('Ploča nakon deformacije')
                
                plt.grid(True)
                plt.savefig('prvi/static/prvi/xyz_plot.png')

                
                Stress[:, ie] = sigmaE.flatten()
                             
                
                sigma_Mises[ie]= np.sqrt(sigmaE[0]**2 + sigmaE[1]**2 - sigmaE[0]*sigmaE[1] + 3*sigmaE[2]**2)
            varijable={
                'XYZ':XYZ,
                'KG':KG,
                'FORCE':FORCE,
                'DISP':DISP,
                'BC':BC,
                'ELDATA':ELDATA,
                'sigmaE':sigmaE,
                'sigma_Mises':sigma_Mises,
                'xyz_plot_path': '/prvi/static/prvi/xyz_plot.png',
                'Stress':Stress,
                
                
            }
            return render(request, 'prvi/LS_rezultati.html', varijable) 

def Preuzimanje(request):
    
    matrica_DISP = np.array2string(DISP, separator=', ')
    matrica_KG = np.array2string(KG, separator=', ')
    matrica_Stress = np.array2string(Stress, separator=', ')

    sadrzaj = f"Matrica pomaka:\n{matrica_DISP}\n\n"
    sadrzaj += f"Matrica krutosti:\n{matrica_KG}\n\n"
    sadrzaj += f"Matrica naprezanja:\n{matrica_Stress}"

    response = HttpResponse(sadrzaj, content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename="Rezultati.txt"'

    return response