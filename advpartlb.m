function[XN,YN] = advpartlb(X0,Y0,U,V,X1,Y1,N,DT)

% ADVPARTLB MATLAB FUNCTION to simulate particles

% uses a Runge Kutta order 2 where:
% X0 is the vector or matrix of longitudes (in degrees)
% Y0 is the vector or matrix of latitudes (in degrees)
% U  is a matrix with the U component of the velocity
% V  is a matrix with the V component of the velocity
% X1 is a vector or matrix with the longitudes of the particles
% Y1 is a vector or matrix with the latitudes of the particles
% N  is the number of segments
% DT is the time increment in hours
%
% If the particle leaves the domain it stays on the border
%
% Made by Jorge Zavala-Hidalgo


DeltaT = DT*3600;
R=6371e+03;                 % Earth's mean radius (Gill)
RH = X0(2)-X0(1);
DeltaTeta=(2*pi/360)*RH;


DelY=R*DeltaTeta;

for ii=1:length(Y0)
    iLat = Y0(ii);
    DelX(ii)=R*DeltaTeta*cos((iLat)*pi/180);
end

[Sx Sy]=size(X0);
if Sx==1 || Sy==1
    [X,Y] = meshgrid(X0,Y0);
else
    X=X0;
    Y=Y0;
end
[Sx Sy]=size(X);

Sx1 =size(X1);

XX1=X1;
YY1=Y1;

XN=XX1;
YN=YY1;

DeltaX = X(1,2) - X(1,1);  % assumes uniform matrix
DeltaY = Y(2,1) - Y(1,1);


[Sx1 Sy1]=size(XX1);
for ii=1:Sx1   % perform the calculations for each vector
    for jj=1:Sy1
        MX=[];
        MY=[];
        %ii
        %jj
        MX(1) = XX1(ii,jj); % coordenadas de las partículas
        MY(1) = YY1(ii,jj);
        k=1;
        HayDatos=0;
        for g=1:N       % calcula para cada segmento del vector
            %g
            J = 1+floor((MX(k) - X(1,1))/DeltaX);
            %I = 1+floor((MY(k) - Y(1,1))/DeltaY);
         IJK = find(Y0 < MY(k));
         I = length(IJK);

            if I >= 1 && I < Sx && J >=1 && J < Sy
                if isnan(U(I,J))~=1 && isnan(U(I+1,J))~=1 ...
                        && isnan(U(I,J+1))~=1 && isnan(U(I+1,J+1))~=1
                    HayDatos=1;
                    %pause
                    IncX = (MX(k) - X(I,J))/DeltaX;
                    IncY = (MY(k) - Y(I,J))/DeltaY;
                    Ua = U(I,J) + (U(I,J+1) - U(I,J)) * IncX;
                    Ub = U(I+1,J) + (U(I+1,J+1) - U(I+1,J)) * IncX;
                    MU(k) = Ua + (Ub - Ua) * IncY;
                    Va = V(I,J) + (V(I,J+1) - V(I,J)) * IncX;
                    Vb = V(I+1,J) + (V(I+1,J+1) - V(I+1,J)) * IncX;
                    MV(k) = Va + (Vb - Va) * IncY;
                    %MX(k+1) = MX(k) + MU(k)*DeltaT;
                    %MY(k+1) = MY(k) + MV(k)*DeltaT;
                    MXInt = MX(k) + MU(k)*DeltaT*cos(MY(k)*pi/180)*(180/(pi*R))/2;  %/(2*R*DelX(I)); % se ubica a delta/2
                    MYInt = MY(k) + MV(k)*DeltaT*(180/(pi*R))/2;  %/(2*R*DelY); % se ubica a delta/2
                    JInt = 1+floor((MXInt - X(1,1))/DeltaX);
                    IInt = 1+floor((MYInt - Y(1,1))/DeltaY);
                    if IInt >= 1 && IInt < Sx && JInt >=1 && JInt < Sy
                        if isnan(U(IInt,JInt))~=1 && isnan(U(IInt+1,JInt))~=1 ...
                                && isnan(U(IInt,JInt+1))~=1 && isnan(U(IInt+1,JInt+1))~=1
                            IncX = (MXInt - X(I,J))/DeltaX;
                            IncY = (MYInt - Y(I,J))/DeltaY;
                            Ua = U(I,J) + (U(I,J+1) - U(I,J)) * IncX;
                            Ub = U(I+1,J) + (U(I+1,J+1) - U(I+1,J)) * IncX;
                            MU(k) = Ua + (Ub - Ua) * IncY;
                            Va = V(I,J) + (V(I,J+1) - V(I,J)) * IncX;
                            Vb = V(I+1,J) + (V(I+1,J+1) - V(I+1,J)) * IncX;
                            MV(k) = Va + (Vb - Va) * IncY;
                            MX(k+1) = MX(k) + MU(k)*DeltaT*cos(MY(k)*pi/180)*(180/(pi*R));   %*DelX(I) DelX(ii)=R*DeltaTeta*cos((iLat)*pi/180);
                            MY(k+1) = MY(k) + MV(k)*DeltaT*(180/(pi*R));  %*DelY
                            k=k+1;
                        end
                    else
                        MX(k+1) = NaN;
                        MY(k+1) = NaN;
                    end
                end
            end
        end
        %      plot(MX,MY,'k');
        %      MX
        %       plot(MX(k+1),MY(k+1),'o','LineWidth',GruesoLinea);
        %ya=1
        %pause
        %      hold on
        if k >= N && HayDatos==1
            %        u(1)=MX(end-1);
            %        v(1)=MY(end-1);
            %        u(2)=MX(end);
            %        v(2)=MY(end);
            %        w=u+v*j;
            %        teta2=atan2(v(2)-v(1),u(2)-u(1));
            %        w2=TPunta*exp(j*(teta2+alfa));
            %        w3=TPunta*exp(j*(teta2-alfa));
            %        w2b(1)=w(2);
            %        w3b(1)=w(2);
            %        w2b(2)=w(2)-w2;
            %        w3b(2)=w(2)-w3;
            %        %plot(w);
            %        plot(w2b,'color',VecColor,'LineWidth',GruesoLinea);
            %        plot(w3b,'color',VecColor,'LineWidth',GruesoLinea) %axis equal
            %        XN(ii)=MX(k);
            %        YN(ii)=MY(k);
            %       clear w w2 w3 w2b w3b
        end
        %ya=MX(k)
        XN(jj)=MX(k);
        YN(jj)=MY(k);
    end
end
%XN
