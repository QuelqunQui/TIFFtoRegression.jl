module TIFF2Regression

############################################################################
# Reading the values of the TIFF according to TIFF v0.6

function FieldValues(ValueOffset::Int64,FileNumber::Array{Int16},FieldCount::Int64,FieldType=4::Int64) #v2
# =========================== Value or Pointer ? ===============================
ValuePointer=0::Int64
if FieldCount==1 # Simplified criterion of "can hold that much info on 4 bytes"
  Value=ValueOffset
else
  ValuePointer=ValueOffset+1 # correction due to the first being 0 and not 1
end

if ValuePointer==0
  DesiredV=Value
else
  #println("$ValuePointer, Count=$FieldCount, Type=$FieldType")
  DesiredV=zeros(1,FieldCount)
  for i=1:FieldCount
    Vi=0
    if FieldType == 2 || FieldType==1 #ASCII or BYTE
      FieldTypeByte=1
    elseif FieldType==3 # Short
      FieldTypeByte=2
    elseif FieldType==4 # Long
      FieldTypeByte=4
    elseif FieldType==5 # Rational
      FieldTypeByte=8
    end

    for jj=1:FieldTypeByte
      Vjj=FileNumber[ValuePointer+FieldTypeByte*(i-1)+jj-1]*256^(jj-1)
      Vi=Vi+Vjj
    end
    DesiredV[i]=convert(Int64,Vi)
  end
end

return DesiredV
end
#############################################################################

# Reading a TIFF (also able to read special Dantec TIFF that are pseudo Std)

function AquiTiff(FilePath::String)
  # Opening File
  FileString=read(FilePath)
  L=length(FileString)
  FileNumber=zeros(Int16,L,1)
  for i=1:L
    NumberI=convert(Int16,FileString[i])
    FileNumber[i]=NumberI
  end
  # ============================ Reading Header ================================
  IFD=FileNumber[5]+1+FileNumber[6]*256+FileNumber[7]*256^2+FileNumber[8]*256^3 # indicates where the Image File Directory starts (i e where the meta data start)
  # +1 in the indexes is a consequence of the starting count from 0 to 1 in tiff ==> JuliaVector
  if FileNumber[1]==73 & FileNumber[2]==73 # Small endian (determines the order of the bites from bigger to smaller)
  # ======================== Reading Metadata: ================================
  # ------------------------- Opening Values ----------------------------------
    Entries=FileNumber[IFD]+256*FileNumber[IFD+1] # Number of Entries in the IFD
    # Opening Values that might be in the Metadata
    Height=0::Int64
    Width=0::Int64
    BPS=0::Int64 #Bit per Sample
    CompType=0::Int64 #Type of Compression
    BlackWhite=0::Int64 #Photometric Interpretation # 0=White is Zero; 1<=Black is Zero
    StripOffsets=0::Int64 #Places where data is stored
    NFragData=0::Int64 # Number of those Places
    Orientation=0::Int64 # 1=Start Top-left
    SPPx=0::Int64 # Sample per Pixel
    RowsPStrip=0::Int64 # Rows per data strip: ceilling(Height/RowsPStrip)=Number of strip (p19)
    StripByteCount=0::Int64 #For each strip the number of bytes in that strip AFTER compression
    XRes=0::Int64 # X-Resolution : nbr of px per Resolution unit (cm/inch see 296)
    YRes=0::Int64 # Y-Resolution
    PlanConfig=0::Int64 # Planar Configuration ???
    ResUnit=0::Int64 # Resolution Unit : abstract (=0), cm (=1) or inch (=2)
    Software=0::Int64 # Software used to create the Tiff
    SBCType=0::Int64
    StripOffsets=[]
    NBytesPPx=0::Int64 # Spec Dantec Dynamics
    # --------------------------------------------------------------------------

    for i=1:Entries # All MetaData there are
      # ----------------------- Initiating Entry--------------------------------
      # Delimiting entry
      StartE=IFD+12*(i-1)+2
      EndE=StartE+11
      EntryN=FileNumber[StartE:EndE] #Taking one entry of 12 hexa bytes
      # Reading Entry
      FieldID=EntryN[1]+EntryN[2]*256
      FieldType=EntryN[3]+EntryN[4]*256
      FieldCount=EntryN[5]+EntryN[6]*256+EntryN[7]*65536+EntryN[8]*16777216
      ValueOffset=EntryN[9]+EntryN[10]*256+EntryN[11]*65536+EntryN[12]*16777216
      # Is either the value of the entry or a pointer towards the place where the values are kept if the entry can't be held on 8bits
      # ----------------------- Switch on FieldID ------------------------------
      if FieldID==254 # NewSubFileType
        NewSubFileType=FieldValues(ValueOffset,FileNumber,FieldCount)
        # 2 : Multipage
        # 4 : Defines a transparency mask
      elseif FieldID==256 # Image Width= n° of column in the Image
        Width=FieldValues(ValueOffset,FileNumber,FieldCount)
      elseif FieldID==257 # Height
        Height=FieldValues(ValueOffset,FileNumber,FieldCount)
      elseif FieldID==258 # Bits/sample ==> number of shade of gray allowable
        BPS=FieldValues(ValueOffset,FileNumber,FieldCount)
      elseif FieldID==259 # Compression Type
        CompType=FieldValues(ValueOffset,FileNumber,FieldCount)
      elseif FieldID==262 # Photometric Interpretation
        # 0=White is Zero; 1<=Black is Zero
        BlackXhite=FieldValues(ValueOffset,FileNumber,FieldCount)
      elseif FieldID==273 # StripOffsets= Position of all data StripOffsets
        NFragData=FieldCount
        StripOffsets=FieldValues(ValueOffset,FileNumber,FieldCount)
      elseif FieldID==274 # Orientation of the Picture
        # 1=Start Top-left
        Orientation=FieldValues(ValueOffset,FileNumber,FieldCount)
      elseif FieldID==277 # Sample/px (3 or more ==> Number of colors allowable)
        SPPx=FieldValues(ValueOffset,FileNumber,FieldCount)
      elseif FieldID==278 # Rows per Strip
        RowsPStrip=FieldValues(ValueOffset,FileNumber,FieldCount)
      elseif FieldID==279 # StripByteCount
        StripByteCount=FieldValues(ValueOffset,FileNumber,FieldCount)
        SBCType=FieldType
      elseif FieldID==282 # X-Resolution
        XRes=FieldValues(ValueOffset,FileNumber,FieldCount)
      elseif FieldID==283 # Y-Resolution
        YRes=FieldValues(ValueOffset,FileNumber,FieldCount)
      elseif FieldID==284 # Planar Configuration
        PlanConfig=FieldValues(ValueOffset,FileNumber,FieldCount)
      elseif FieldID==285
        PageName=FieldValues(ValueOffset,FileNumber,FieldCount,FieldType)
      elseif FieldID==296 # Resolution Unit (cm/inch)
        ResUnit=FieldValues(ValueOffset,FileNumber,FieldCount)
      elseif FieldID==297
        PageNumber=FieldValues(ValueOffset,FileNumber,FieldCount,FieldType)
      elseif FieldID==305 # Software
        #Software=FieldValues(ValueOffset,FileNumber,FieldCount)
      elseif FieldID==32997
        #Specific to DynamicStudio, Dantec
        # Number of bytes/ px? TBC
        NBytesPPx=FieldValues(ValueOffset,FileNumber,FieldCount)
        if SPPx==0
          SPPx=2 #NBytesPPx/8 # Works?
          SPPx=convert(Int64,SPPx)
        end
      end  # Metadata Aquisition
    end
    EndIFD=IFD+12*(Entries-1)+13
  elseif FileNumber[1]==77 & FileNumber[2]==77 #big endian
    println("Big Endian Still to be implemented, nearly same as small endian")
  end
  # ======================= Creating the Img Matrix ============================
  # ----------------------- Finding and Ordering Px ---------------------------
  PixelString=[]
  PixelNumber=ones(1,length(FileNumber))
  CountLeng=0
  for i=1:NFragData
    PxStrip=FileNumber[convert(Int64,StripOffsets[i]):convert(Int64,StripOffsets[i]+StripByteCount[i])]
    l=length(PxStrip)
    PixelNumber[1+1+CountLeng:1+l+CountLeng]=PxStrip[1:l]
    CountLeng=CountLeng+l
  end
  PixelNumber=PixelNumber[1:CountLeng] # cutting the exedentary ones (due to the last strip being sometimes shorter than the others)
  # ------------------------ Making Order in the Px ----------------------------
  ImgArray=zeros(Int64,Height,Width,SPPx)
  for i=1:Height
    for  j=1:Width
      for k=1:SPPx
        if NBytesPPx==0
          Pxijk=((i-1)*Width+(j-1))*SPPx+k
        else # Why ? 171 dtmed empirically for Dynamic Studio Tif)
          #Pxijk=((i-1)*Width+(j-1))*SPPx+k
          Pxijk=((i-1)*Width+(j+convert(Int64,floor(171*i/Height))-1))*SPPx+k # ???????????????????????????????????????????????????
        end
        ImgArray[i,j,k]=PixelNumber[Pxijk]
      end
    end
  end
  #--------------------- Making it Gray Matrix ---------------------------------
  if NBytesPPx==0 # Suposedly "Normal Img or Color Img"
    #GrayArray=sum(ImgArray,3)/SPPx
    GrayArray=ImgArray[:,:,3]
  else
    GrayArray=ImgArray[:,:,1]*256+ImgArray[:,:,2] # Dantec Spec
  end
return GrayArray[:,:,1]
end

###########################################################################

# Loading a Serie of TIFF from Dantec Dynamics

function AquiSerieD(Folder::String, File0::String, NbrImg=189::Int64)
  # ========================== Making Img Path =================================
  # ----------------Cutting File0 (Dantec DynamicsStudio Naming)----------------
  j=1
  while File0[j]!='.'
    j=j+1
  end
  p1=j
  #println(p1)
  BaseName=File0[1:p1-1] # Input by User exporting Img
  j=j+1
  while File0[j]!='.'
    j=j+1
  end
  p2=j
  FileSerie=File0[p1+1:p2-1] #Auto Serie Numbering
  NumberInSerie=File0[p2+1:end] # Auto Img Numbering
  N=NbrImg-1 # Numbering Starts at 0 not 1
  ImgArray=[]
  OneAdress=[]
  # ------------------------ Making Image Number String -------------------------
  for i=0:N
    # Finding k
    LNIS=length(NumberInSerie)
    k=LNIS
    while floor(i/(10.0^k))<1 && k>=0
      k=k-1
    end
    # completing the NumberIn Serie String
    NIS=""
    ik=i
    for kk=LNIS:-1:1
      if kk>k+1
        NIS="$(NIS)0"
      else
        V=convert(Int64,floor(ik/(10^(kk-1))))
        NIS="$(NIS)$V"
        ik=ik-V*10^(kk-1)
      end
    end
    OneAdress="$(Folder)\\$(BaseName).$(FileSerie).$(NIS).tif"
    #println(OneAdress)
    # ========================= Aquiring Images ==================================
    OneImg=AquiTiff(OneAdress)
    if isempty(ImgArray)
      ImgArray=zeros(Int64,size(OneImg,1),size(OneImg,2),N+1)
      ImgArray[:,:,1]=OneImg
    else
      ImgArray[:,:,i+1]=OneImg
    end
  end
return ImgArray # to assure Julia wont print the Array
end
#############################################################################

# Loading a Serie of TIFF automaticaly renamed with Windows 10

function AquiSerieW(Folder::String, File0::String, NbrImg=189::Int64)
  # ========================== Making Img Path =================================
  # ----------------Cutting File0 (Dantec DynamicsStudio Naming)----------------
  j=1
  while File0[j]!='('
    j=j+1
  end
  p1=j
  #println(p1)
  BaseName=File0[1:p1-1] # Input by User + 1 space (std windows renumerotation)
  j=j+1

  ImgArray=[]
  OneAdress=[]
  # ------------------------ Making Image Number String -------------------------
  for i=1:NbrImg
    OneAdress="$(Folder)\\$(BaseName)($i).tif"
    #println(OneAdress)
    # ========================= Aquiring Images ==================================
    OneImg=AquiTiff(OneAdress)
    if isempty(ImgArray)
      ImgArray=zeros(Int64,size(OneImg,1),size(OneImg,2),NbrImg)
      ImgArray[:,:,1]=OneImg
    else
      ImgArray[:,:,i]=OneImg
    end
  end
return ImgArray # to assure Julia wont print the Array
end

#############################################################################

# Creating a Special type to stock the start of a line and its values
type Line # Creates type Line
  YValues::Vector
  XStart::Int64
end
##############################################################################

# Detecting the Lines in an image matrix (also retourning a Binary matrix)

function DetectLine(Img::Matrix)
  # only for bicolor ImgLine
  BiMat=zeros(Bool,size(Img))
  BiMat0=zeros(Bool,size(Img))
  Max=maximum(Img)
  min=minimum(Img)
  moy=(Max+min)/2
# --------------------------BiColoring----------------------------------
  k=0
  for i=1:size(Img,1) #BiLeveling
    for j=1:size(Img,2)
      if Img[i,j]<=moy
        BiMat[i,j]=0
        BiMat0[i,j]=0
      else
        BiMat[i,j]=1
        BiMat0[i,j]=1
      end
    end
  end

SumMat=sum(BiMat)
H=size(BiMat,1)
L=size(BiMat,2)
NPx=H*L
if SumMat<NPx/2
  # 1=Line, 0=Background
  LinePx=true
  BckGrnd=false #BckGrnd LinePx
else
  # 0=Line, 1=Background
  LinePx=false
  BckGrnd=true
end
CountLinePx=abs(SumMat-NPx) # nbr of Px that are part of a line
Lines=[] # need Array to stock the XStart
Xcount=1::Int64
Split=false
Ymean2=H
YmeanM1=H
Xcount2=1
while Xcount<L+1
  x=Xcount
 # 1 Line code  !
  if Split
    y=Ymean2
    x=Xcount2
    Xcount=Xcount-1
    Split=false
  else
    y=H
    x=Xcount
  end
  ConLine=true
  # ========================= Following the Line =========================
  while y>=2 && ConLine# going up to top
    # find [Background; Line]
    YmeanV=[]
    if [LinePx; BckGrnd]==[BiMat[y-1,x]; BiMat[y,x]] # Line Start Detect
      #Xcount=Xcount-1 #Stay in the same column and check further
      # Special Start in case of Split
      if YmeanM1 != H
        push!(YmeanV,YmeanM1)
        XStart=x-1
        YmeanM1=H
      else #End Spec Split
        XStart=x
      end
      y=y-1
      YStart=y
      # Start of a Line
      while BiMat[y,x]==LinePx
        y=y-1
      end
      YEnd=y+1
      # following that one line
      Width=abs(YStart-YEnd)+1
      BiMat[YEnd:YStart,x]=repmat([BckGrnd],Width,1)
      Ymean=(YStart+YEnd)/2
      push!(YmeanV,Ymean)
      while ConLine && x<L
        x=x+1
        if YStart<H && YEnd>1
          if sum(BiMat[YEnd-1:YStart+1,x].==repmat([BckGrnd],Width+2,1))==Width+2
            # Line Stop
            ConLine=false
            YStart=0
            YEnd=0
            x=Xcount
            Xcount=Xcount-1
          else
            # Line continue
            if BiMat[YEnd-1:YStart+1,x]==[BckGrnd;repmat([LinePx],Width,1);BckGrnd]
             # Line constant in Value and Width
            else
              # Line continue but changed
              if (BiMat[YStart:YStart+1,x]==[LinePx;LinePx])||(BiMat[YEnd-1:YEnd,x]==[LinePx;LinePx])
                # Line increase or decrease
                if (BiMat[YStart:YStart+1,x]==[LinePx;LinePx]) && !(BiMat[YEnd-1:YEnd,x]==[LinePx;LinePx])
                  # Decreasing
                  YEnd=YEnd+minimum(findin(BiMat[YEnd:YStart,x].==repmat([LinePx],Width,1),true))-1
                  yy=YStart
                  while BiMat[yy+1,x]==LinePx && yy<H
                    yy=yy+1
                  end
                  YStart=yy
                elseif !(BiMat[YStart:YStart+1,x]==[LinePx;LinePx]) && (BiMat[YEnd-1:YEnd,x]==[LinePx;LinePx])# Increasing
                  YStart=YEnd+maximum(findin(BiMat[YEnd:YStart,x].==repmat([LinePx],Width,1),true))-1
                  # YEnd
                  yy=YEnd
                  while BiMat[yy-1,x]==LinePx && yy>1
                    yy=yy-1
                  end
                  YEnd=yy
                elseif (BiMat[YStart:YStart+1,x]==[LinePx;LinePx]) && (BiMat[YEnd-1:YEnd,x]==[LinePx;LinePx]) # Increse in up and down
                  if sum(BiMat[YEnd-1:YStart+1,x].==LinePx)==(YStart+1-(YEnd-1)+1)# No hole
                    # Increse in Width in both side
                    ydown=YStart
                    while BiMat[ydown+1,x]==LinePx && ydown<H
                      ydown=ydown+1
                    end
                    YStart=ydown
                    yup=YEnd
                    while BiMat[yup-1,x]==LinePx && yup>1
                      yup=yup-1
                    end
                    YEnd=yup
                  else # there is a hole
                    # Split in two curves
                    Split= true
                    PosHole=findin(BiMat[YEnd-1:YStart+1,x].==BckGrnd,true)
                    # Down Part Continues
                    YEnd=(YEnd-1+minimum(PosHole))-1 # (below hole)
                    ydown=YStart
                    while BiMat[ydown+1,x]==LinePx && ydown<H
                      ydown=ydown+1
                    end
                    YStart=ydown
                    # Up Going part
                    YStart2=(YEnd-1+maximum(PosHole))+1
                    Ymean2=YStart2+1
                    Xcount2=x
                    YmeanM1=Ymean # kept for next line
                  end
                end
              else
                #Line Width decrease within previous limits
                PxPos=findin(BiMat[YEnd:YStart,x].==repmat([LinePx],Width,1),true)
                minPxPos=minimum(PxPos)
                MaxPxPos=maximum(PxPos)
                YStart=YEnd+MaxPxPos-1
                YEnd=YEnd+minPxPos-1
              end
            end
          end
        else
          # border effect
        end
        if !ConLine
          Width=0
        else
          Width=abs(YStart-YEnd)+1
          BiMat[YEnd:YStart,x]=repmat([BckGrnd],Width,1)
          Ymean=(YEnd+YStart)/2
          push!(YmeanV,Ymean)
        end
      end #start of a line detected
    else
      y=y-1 # No Line, go on climbing
    end
    if length(YmeanV)>1 #!isempty(YmeanV)
      Line1=Line(H+1-YmeanV,XStart) # Assuming left bottom axis =0
      #then complete Lines
      push!(Lines,Line1)
      ConLine=false # end current line
    end
  end
    Xcount=Xcount+1
end
# from bottom to top
# 1 First column -> detect LinePx
# 1b if no LinePx go to the right
# 2 continue climbing until Px!=LinePx
# 3 from there moove to the right following the line
# 3a check same Px and Px+1 and Px-1
# 3b if Px+-1 check also Px+-2, etc
#    if no Px: line stop
# 3c repeat 3 until line stop
# 4 remoove selected Px
# 5 repeat operation until top right
return BiMat0, Lines
end

#############################################################################

# Allowing to plot the new type created, eventually on the matrix (thanks to PyPlot)

function PlotLine(Lines::Array{Any},H=0::Int64)
  for i=1:length(Lines)
    XS=Lines[i].XStart
    YV=Lines[i].YValues
    if H==0
      plot(XS:XS+length(YV)-1,YV,"-")
    else
      plot(XS-1:XS-1+length(YV)-1,H-YV,"-")
    end
  end
end

############################################################################

# Creating the next function for the Multiple Regression Function below

function NextFun(Xin::Vector,RegFuns="Poly"::String,i=1::Int64)
  n=length(Xin)
  if RegFuns=="Poly"
    Xout=Xin.^(i)
    XT=("Poly Order $i")
  elseif RegFuns=="Trigo"
    if isinteger(i/2)
      Xout=cos.(Xin/maximum(Xin)*(i/2))
      XT=("Standirdize cos($i x/2)")
    else
      Xout=cos.(Xin/maximum(Xin)/((i+1)/2))
      XT=("Standirdize cos($(1/(i+1)) x/2)")
    end
  elseif RegFuns=="Frac"
    F=findin(Xin,0)
    for j=1:length(F)
      X[F[j]]=mean(abs(Xin))*10^-5
    end
    Xout=1./Xin.^i
    if abs(minimum(Xout))< 2e-15 # all becoming too small
      println("Elements becoming too small X=ones")
      Xout=ones(size(Xin))
      XT=("Ones")
    else
      XT=("Frac Order $i")
    end
  elseif RegFuns=="Expo"
    Xout=exp.(i*Xin/maximum(Xin))
    XT=("Standardize exp($i x)")
  elseif RegFuns=="Loga"
    Xout=log.(abs(Xin+(i-2)*sqrt(maximum(Xin))))
    XT=("log(x+$(i-2) sqrt(Max))")
  elseif RegFuns=="Mixte"
    NFam=5 #nbr of family used (code written could be improoved)
    if isinteger((i+NFam-1)/NFam)
      k=(i+NFam-1)/NFam
      Xout=Xin.^k
      XT=("Poly Order $k") #Poly
    elseif isinteger((i+NFam-2)/NFam)
      k=(i+NFam-2)/NFam
      F=findin(Xin,0)
      for j=1:length(F)
        X[F[j]]=mean(abs(Xin))*10^-5
      end
      Xout=1./Xin.^k
      if abs(minimum(Xout))< 2e-15 # all becoming too small
        println("Elements becoming too small X=rand")
        Xout=rand(size(Xin))
      else
        XT=("Frac Order $k")
      end #Frac
    elseif isinteger((i+NFam-3)/NFam)
      k=(i+NFam-3)/NFam
      #NTrigo=4
      if isinteger(k/2)
        Xout=cos.(Xin/maximum(Xin)*(k/2))
        XT=("Standirdize cos($k x/2)")
      else
        Xout=cos.(Xin/maximum(Xin)/(k+1)/2)
        XT=("Standirdize cos($(1/(k+1)) x/2)")
      end #Trigo
    elseif isinteger((i+NFam-4)/NFam)
      k=(i+NFam-3)/NFam
      Xout=exp.(k*Xin/maximum(Xin))
      XT=("Standardize exp($k x)") #Expo
    elseif isinteger((i+NFam-5)/NFam)
      k=(i+NFam-5)/NFam
      Xout=log.(abs(Xin+(k-2)*sqrt(maximum(Xin))))
      XT=("log(x+$(k-2) sqrt(Max))")
    end
  elseif RegFuns=="Custom" # Would their be another interesting family
    println("What next fun do you want ? (Existing type: Poly, Trigo, Frac, Expo or Loga) \n Your choice is :")
    SwCu=input() #Switch Custom
    if SwCu=="Poly"
      println("Order of the wanted term : y=x^{Order} ? \n Order=")
      u=parse(Float64,input())
      if isinteger(u)
        Xu=Xin
      else
        Xu=Xin
        Xu[find(Xu.<0)]=zeros(size(find(Xu.<0)))
        println("All smaller than 0 elements have been replaced by zero for Poly Order $u")
      end
      Xout=Xu.^(u)
      XT=("Poly Order $u")
    elseif SwCu=="Trigo"
      println("Normal or Inverse ? \n Your choice is :")
      TTrig=input()
      if TTrig=="Normal"
        println("What type do you want ? (Possible choices: cos, sin, tan or cot) \n Your choice is :")
        TrigN=input()
        if TrigN=="sin"
          pritnln("What Angulat fequency (omega) do you want ? \n omega =")
          omega=parse(Float64,input())
          println("What Dephasage (phi) do you want ? \n phi =")
          phi=parse(Float64,input())
          Xout=sin.(omega*Xin+phi)
          XT=("sin($omega x+$phi)")
        elseif TrigN=="cos"
          pritnln("What Angulat fequency (omega) do you want ? \n omega =")
          omega=parse(Float64,input())
          println("What Dephasage (phi) do you want ? \n phi =")
          phi=parse(Float64,input())
          Xout=cos.(omega*Xin+phi)
          XT=("cos($omega x+$phi)")
        elseif TrigN=="tan"
          pritnln("What Angulat fequency (omega) do you want ? \n omega =")
          omega=parse(Float64,input())
          println("What Dephasage (phi) do you want ? \n phi =")
          phi=parse(Float64,input())
          Xout=tan.(omega*Xin+phi)
          XT=("tan($omega x+$phi)")
        elseif TrigN=="cot"
          pritnln("What Angulat fequency (omega) do you want ? \n omega =")
          omega=parse(Float64,input())
          println("What Dephasage (phi) do you want ? \n phi =")
          phi=parse(Float64,input())
          Xout=cot.(omega*Xin+phi)
          XT=("cot($omega x+$phi)")
        end
      elseif TTrig=="Inverse"
        println("What type do you want ? (Possible choices: acos, asin, atan or acot) \n Your choice is :")
        TrigN=input()
        if TrigN=="asin"
          pritnln("Vector is standardize (Divided by Max)")
          println("f ? : asin(XStd/f) (abs(f)>1)? \n f=")
          f=parse(Float64,input())
          if abs(f)<1
            f=1
          end
          Xout=asin.(Xin/maximum(abs(Xin))/f)
          XT=("asin(XStd/$f)")
        elseif TrigN=="acos"
          pritnln("Vector is standardize (Divided by Max)")
          println("f ? : acos(XStd/f) (abs(f)>1)? \n f=")
          f=parse(Float64,input())
          if abs(f)<1
            f=1
          end
          Xout=acos.(Xin/maximum(abs(Xin))/f)
          XT=("acos(XStd/$f)")
        elseif TrigN=="atan"
          println("What Factor f ? (atan(fX+t)) \n f=")
          f=parse(Float64,input())
          println("What terme ? \n t=")
          t=parse(Float64,input())
          Xout=atan.(f*Xin+t)
          XT=("atan(fX+t)")
        elseif TrigN=="acot"
          println("What Factor f ? (acot(fX+t)) \n f=")
          f=parse(Float64,input())
          println("What terme ? \n t=")
          t=parse(Float64,input())
          Xout=acot.(f*Xin+t)
          XT=("acot(fX+t)")
        end
      end
    elseif SwCu=="Frac"
      println("Order of the wanted term : y=x^-{Order} ? \n Order=")
      u=parse(Float64,input())
      if isinteger(u)
        Xu=Xin
      else
        Xu=Xin
        Xu[find(Xu.<=0)]=ones(size(find(Xu.<0)))
        println("All smaller than or equal to 0 elements have been replaced by one for Frac Order $u")
      end
      Xout=1./(Xu.^(u))
      XT=("Frac Order $u")
    elseif SwCu=="Expo"
      println("Order of the wanted Exponential ? : y=exp(x*{Order}) ? \n Order=")
      u=parse(Float64,input())
      Xout=exp.(u*Xin)
      XT=("exp($u x)")
    elseif SwCu=="Loga"
      println("Terme of the wanted Logarithme ? : y=log|x+t|) ? \n t=")
      t=parse(Float64,input())
      Xout=log.(abs(Xin+t))
      XT=("log(x+t)")
    end
  end
return Xout,XT
end

###############################################################################

# Statistics discrabing the Regression using ANOVA analysis thanks to distributions
# functions

function StatReg(Xin,Yin,XMat,XName,p)
  Yb=mean(Yin)
  SquareX=(transpose(XMat)*XMat)
  n=length(Yin)
  Bh=(SquareX^-1)*transpose(XMat)*Yin # Estimator for Inv(Mat) ?
  Yh=XMat*Bh
  SSE=sum((Yin-Yh).^2)
  SSR=sum((Yh-Yb).^2)
  SSTot=sum((Yin-Yb).^2)
  R2=SSR/SSTot # the biggest the best
  dfR=p-1
  MSR=SSR/dfR
  dfE=n-p
  MSE=SSE/dfE
  F=MSE/MSR
  # ---------------------- Calc p-Value Model --------------------------
  FD=FDist(dfR,dfE)
  pModel=cdf(FD,F) #Smallest the best
  # p- Values indie
  TD=TDist(dfE)
  pVect=ones(size(Bh))
  for i=1:length(Bh)
    Yh=XMat[:,i]*Bh[i]
    SSE=sum((Yin-Yh).^2)
    if XName[i]=="Constant"
      SX2=sum((Xin-mean(Xin)).^2)
      SEBh=sqrt(SSE*(1/n+mean(Xin)^2/SX2))
    else
      X=XMat[:,i]
      SX2=sum((X-mean(X)).^2)
      SEBh=sqrt(SSE/SX2)
    end
    t=Bh[i]*sqrt(dfE)/SEBh
    pVar=2*(1-cdf(TD,t))
    pVect[i]=pVar
  end

  return pVect,pModel,Bh
end

#############################################################################

# Main function of Multiple Regression implementing ANOVA tests
# and a user friendly approach to regression
############################################################################
end # module
