classdef Geometric_3D_Surface_Reconstruction_from_1D_Laser_Data_v1_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        CorningIncLaserGaugeDataAnalysisUIFigure  matlab.ui.Figure
        AngleNo                    matlab.ui.control.NumericEditField % No. of angles for 3D cylindrical display
        AngleLabel                 matlab.ui.control.Label % angle
        Label2                     matlab.ui.control.Label % ui label
        ExportresultsButton        matlab.ui.control.Button
        Label                      matlab.ui.control.Label % cylindrial height
        PlaneheightDropDown        matlab.ui.control.DropDown
        PlaneheightDropDownLabel   matlab.ui.control.Label
        PartName                   matlab.ui.control.DropDown % part
        PartIDLabel                matlab.ui.control.Label % part ID input
        PartRadius                 matlab.ui.control.NumericEditField % radius of the cylindrical part
        PartHeight                 matlab.ui.control.NumericEditField % height of the part
        plot                       matlab.ui.control.Button %function button to plot figure
        ImportdatafileButton       matlab.ui.control.Button %button for import data
        Scalefactor                matlab.ui.control.NumericEditField % scale ratio for display
        ScalefactorEditFieldLabel  matlab.ui.control.Label
        HeightEditFieldLabel       matlab.ui.control.Label
        RadiusEditFieldLabel       matlab.ui.control.Label 
        UIA3d                      matlab.ui.control.UIAxes %3D plot
        UIA2d                      matlab.ui.control.UIAxes %2D plot
        UIA2dcontour               matlab.ui.control.UIAxes %contour plot
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    properties (Access = public)
        % Property description
        VIDECK =[];
        PartID=strings;
        planeheight=[];
        ID; %Part ID for plot
        nangle;
        NP; %number of planes
        Radius;
        Height;
        ScaleF;
        HeightEdit;
        partflag=2;
        heightflag=2;
        LG_Scaled=[];
        spiralpara=[]; %spiral function parameters
        dataflag=2;
        intrp='Natural-Natural'; % interpolation method
        Pax; % Polar axes
    end
    
methods (Access = public)
    function Readfile_ExtractData (app, event) %read in laser gauge measurement data file
            [fname, path] = uigetfile ('*.csv*');
            if (fname)
                app.Label.Text = fname;
                inpuFilePath = fullfile(path,fname);
                prop=readtable(inpuFilePath,'ReadVariableNames',true,'HeaderLines',1);
                app.NP = table2array(prop(1,'ContourPlanes'))
                filecolumn=max(size(prop(1,:)));
                filerow =max(size(prop(:,1)));
                columnpart =find(strcmpi(prop.Properties.VariableNames,'Piece')); %*findColNumber*(prop, 'Piece');
                columnpheight = find(strcmpi(prop.Properties.VariableNames, 'PlaneHeight'));
                app.planeheight= table2array(prop(1:app.NP,columnpheight));
                columnvidek = find(strcmpi(prop.Properties.VariableNames, 'VidekPoints1'));
                app.nangle = app.AngleNo.Value;   %filecolumn-columnvidek+1; 
                ipart=0;
                partID=strings;
                for i=1:app.NP:filerow
                        ipart=ipart+1; %
                        app.PartID (ipart) = table2array (prop(i,columnpart)); %1x50 strings
                        rowj=(ipart-1)*app.NP+1;
                        app.VIDECK(ipart,1:app.NP,1:app.nangle) = table2array(prop(rowj:rowj+app.NP-1,columnvidek:columnvidek+app.nangle-1)); % 50x15x24 double
                        app.VIDECK(ipart,1:app.NP,app.nangle+1) = table2array(prop(rowj:rowj+app.NP-1,columnvidek)); % 50x15x25
                end
                app.PartName.Items = app.PartID;
                txt0 (1)={'All'};
                for i=2:app.NP+1
                    txt0 (i)= {['h',num2str(i),'=',num2str(round(app.planeheight(i-1),2))]};
                end
                app.PlaneheightDropDown.Items=txt0;
                app.HeightEdit=0; % median(app.planeheight);
                app.dataflag=1;
                app.ScaleF = app.Scalefactor.Value;
                app.Height = app.PartHeight.Value;   
                app.Radius = app.PartRadius.Value;   
                
    else
        disp('User selected Cancel')
        return 
    end    
end %end of function

function resultsexport(app,event) % export results out as csv
    LGD=[]; % 
    T =[]; % numerical or cell?? array to store results FormalTable({'Part_ID','Mean','Height','Weight'}); % (Part_ID,Mean,Height,Weight);
    T{1,1}='Part ID';
    T{1,2}='Deviation_Mean';
    T{1,3}='Deviation_Max';
    T{1,4}='Deviation_Min';
    T{1,5}='Deviation_height';
    T{1,6}='Deviation_STD'; % (Part_ID,Mean,Height,Weight)
    h=app.planeheight;
    nr=(app.nangle+1);
    X= [0:2*pi/app.nangle:2*pi]; % angle in radius
    nonodes_tan = 361; % 73, 5 degree per data
    dtheta=double(2*pi/(nonodes_tan-1));  %Radians at each measurement increment
    dL=dtheta*app.Radius;
    nonodes_axial = int16(app.Height/dL)+1;
    deltah = app.Height/(double(nonodes_axial)-1); %strainge that the "double" has to be added otherwise it returend as 0 (int value)
    ntotal=double(nonodes_axial)*double(nonodes_tan);
    nodes = zeros (ntotal,3);
    total_theta=[0:dtheta:360];
    nodes(:,1)=[1:1: ntotal];
    for i = 1 : ntotal
        nodes(i,2)=double((i-floor((double(i)-1)/nonodes_tan)*nonodes_tan-1))*dL; %along unroll clockwise direction
        nodes(i,3)= double(floor((double(i)-1)/nonodes_tan))*deltah; %along height
    end    
    for i = 1 : app.NP*nr
        row=floor((i-1)/nr) + 1;
        column =i-(row-1)*nr;
        MECoords(i,2)=X(column)*app.Radius;% length as X value, unroll clockwise direction %double((i-floor((double(i)-1)/nr)*nr-1))*dr;    
        MECoords(i,3)= h(row);
    end
    partnumber=max(size(app.PartID));
    
    i_mehtod=2 %fit/exptrapolation method, 1: natural/linear, 2: spiral in radial, and pchit in height    
    [xi,yi] = meshgrid(0: 2*pi*app.Radius/double(nl):2*pi*app.Radius,0:app.Height/double(nh):app.Height);
    
    
for ipart=1:partnumber
    LGD (:,:) = app.VIDECK(ipart,:,:);
    for i = 1 : app.NP*nr
        row=floor((i-1)/nr) + 1;
        column =i-(row-1)*nr;
        VIDEK (i) = LGD (row, column);
    end    
    %interpolation-extrapolation
    if (i_mehtod==1) %
        F = scatteredInterpolant(MECoords(:,2:3),VIDEK','pchip','linear'); %VIDEK' must be column wise
        nVIDEK=F(nodes(:,2),nodes(:,3));    
        nl=500;
        nh=int16(nl*app.Height/(2*pi*app.Radius));
        %F.Method='natural';
        F.Method='pchip';
        F.ExtrapolationMethod='linear';    
        zi=F(xi,yi);    % Solve the interpolant of VIDECK at desired grid points using the natural-linear interpolation/extrapolation method
    else
        
    end
           
    Deviation_Mean = mean2(zi); %roughness Raa
    Deviation_STD =std2(zi);
    Deviation_Max = max(max(zi));
    Deviation_Min = min(min(zi));
    Deviation_height = Deviation_Max-Deviation_Min;
    Deviation_Rq = sqrt(mean(zi.^2)); %Root mean squared
    Deviation_Rsk = (1/Deviation_Rq^3)*mean(zi.^3); % skewness
    Deviation_Rku = (1/Deviation_Rq^4)*mean(zi.^4); % skewness
    
    T(ipart+1,1)= cellstr(app.PartID(ipart));
    T(ipart+1,2)=num2cell(Deviation_Mean);
    T(ipart+1,3)= num2cell(Deviation_Max);
    T(ipart+1,4)= num2cell(Deviation_Min);
    T(ipart+1,5)= num2cell(Deviation_height);
    T(ipart+1,6)= num2cell(Deviation_STD);
end % end of loop part ID
    [file,path] = uiputfile('results.xlsx');
    resultFilePath = fullfile(path,file);
    writecell(T, resultFilePath);
    app.Label2.Text ='';
end % end of function - results export



function mainplot (app,event)    
    app.ID = char(app.PartName.Value);    
    app.intrp = 'Natural-Linear'; % app.InterpoextrapolationDropDown.Value;
    LGD=[];
    if length(app.ID)==7 & app.ID(1:7)=='Part ID'
        display ('You must pick up a Part ID for analysis')
        return;
    else
        partnumb = find(app.ID==app.PartID); % compare between double and string
        LGD (:,:) = app.VIDECK(partnumb,:,:);
    end 
    h=app.planeheight;
    nangle=app.nangle;
    NP=app.NP;
    height=app.Height;
    Radius=app.Radius;
    ScaleF=app.ScaleF;
    nr=(nangle+1);
    X= [0:2*pi/nangle:2*pi]; %angle in radius 
    for i = 1 : NP*nr % Np = 15 nr= 25, number of plane, number of radial points (24 + 1)
        row=floor((i-1)/nr) + 1;
        column =i-(row-1)*nr;
        MECoords(i,2)=X(column)*Radius;% length as X value, unroll clockwise direction %double((i-floor((double(i)-1)/nr)*nr-1))*dr;    
        MECoords(i,3)= h(row);
        VIDEK (i) = LGD (row, column);
    end
    nonodes_tan = 361; % 73, 5 degree per data
    dtheta=double(2*pi/(nonodes_tan-1));  %radians at each measurement increment
    dL=dtheta*Radius;
    nonodes_axial = int16(height/dL)+1;
    deltah = height/(double(nonodes_axial)-1); %strainge that the "double" has to be added otherwise it returend as 0 (int value)
    ntotal=double(nonodes_axial)*double(nonodes_tan);
    nodes = zeros (ntotal,3);
    total_theta=[0:dtheta:360];
    nodes(:,1)=[1:1: ntotal];
    for i = 1 : ntotal
        nodes(i,2)=double((i-floor((double(i)-1)/nonodes_tan)*nonodes_tan-1))*dL; %along unroll clockwise direction
        nodes(i,3)= double(floor((double(i)-1)/nonodes_tan))*deltah; %along height
    end
    
    i_interpo=1
    
    nl=500;
    nh=int16(nl*height/(2*pi*Radius));
    [xi,yi] = meshgrid(0: 2*pi*Radius/double(nl):2*pi*Radius,0:height/double(nh):height);
    if i_interpo==1
        %interpolation-extrapolation
        F = scatteredInterpolant(MECoords(:,2:3),VIDEK','natural','linear'); %VIDEK' must be column wise
        nVIDEK=F(nodes(:,2),nodes(:,3));
        F.Method='natural';
        F.ExtrapolationMethod='linear'; % method testing here,'pchip' doesnot work
        zi=F(xi,yi);    % Solve the interpolant of VIDECK at desired grid points    
    else
        for i = 1 : NP % Np = 15 nr= 25, number of plane, number of radial points (24 + 1)
        xx(i,:)=X
        end
        for j=1:nr
        yy(:,j)=h
        end
        zi = interp2(xx,yy,LGD,xi,yi,'makima');
        figure
        surf(xi,yi,zi)
    end
    
    theta=[0:2*pi/(nangle):2*pi];  % 1 row, radians at each measurement increment
    
if app.heightflag==2    %analysis of specific height 'off'
    outplot=2;    
    ndata=NP; %length(LG(:,1));    
    scale=ScaleF;    
    app.LG_Scaled = Radius + LGD (1:NP,:)*scale;  %Scaled Deformation app.LG_Scaled = R + LG.*scale;  %Scaled Deformation   
    % Unrolled Contour surface Plot (2D)
    hold(app.UIA2dcontour,'on');
   [C,h]=contourf(app.UIA2dcontour,int16(xi/Radius*180/pi),yi,zi, 100,'LineStyle','none') % );'c-',,%plot the unrolled 2D contour surface, xi: circular angle, yi: height, zi: deviation

    grid (app.UIA2dcontour,'on');    
    Deviation_Rq = sqrt(mean2(zi.^2)); %Root mean squared
    Deviation_Rsk = (1/Deviation_Rq^3)*mean2(zi.^3); % skewness
    txt1 = {['Mean=' num2str(round(mean2(zi),5))],['Max=' num2str(round(max(max(zi)),5))],['Min=' num2str(round(min(min(zi)),5))],['STD=' num2str(round(std2(zi),5))],['Skewness=' num2str(Deviation_Rsk)]};
    text(app.UIA2dcontour,180,app.Height/2,txt1,'FontSize',10); %,'HorizontalAlignment','left')
    xticks(app.UIA2dcontour,[0 30 60 90 120 150 180 210 240 270 300 330 360]);
    xlim(app.UIA2dcontour,[0 360]);
    set(app.UIA2dcontour,'XDir', 'reverse');
    ylim(app.UIA2dcontour,[0 app.Height]);
    colorbar (app.UIA2dcontour); % title(hcb,'Radius deviation (in)')
    hold(app.UIA2dcontour,'off');
    app.partflag=1;
    if app.partflag   % partflag=1: already have part number;2: no part number for analysis
        theta1=[0:2*pi/360:2*pi];
        hei=app.HeightEdit;
        app.PlaneheightDropDown.Value = "All";
        hold(app.UIA2d,'on');
        plot(app.UIA2d,Radius*cos(-theta1),Radius*sin(-theta1),'LineWidth',1,'Color','r','LineStyle',':') % plot nomial
        i_plane=1;
        while i_plane< app.NP+1
            plot(app.UIA2d,app.LG_Scaled(i_plane,:).*cos(-theta),app.LG_Scaled(i_plane,:).*sin(-theta),'LineWidth',1,'LineStyle','-'); % plot all measurement data 
            i_plane=i_plane+1;
        end
        legendCell =['Nomial';cellstr(num2str(round(app.planeheight,2), 'h=%.2f'))];
        txt = {['Measurement:'],['Mean=' num2str(mean2(LGD))],['Max=' num2str(max(max(LGD)))],['Min=' num2str(min(min(LGD)))],['STD=' num2str(std2(LGD))]};
        grid (app.UIA2d, 'on');
        for i=1:12
            txt0 ={(i-1)*30};
            text(app.UIA2d,0.9*Radius*cos(-(i-1)*pi/6),0.9*Radius*sin(-(i-1)*pi/6),txt0,'Rotation',180-(i-1)*30,'Color','red','FontSize',10);
        end
        legend (app.UIA2d,legendCell,'Location','southeastoutside','color','none','fontSize',6); %'Orientation','horizontal'
        text(app.UIA2d,-app.Radius/2.5,0.1*app.Radius,txt,'FontSize',10); %,'HorizontalAlignment','left')
        hold(app.UIA2d,'off');
    end
    
    % 3d surface Plot
    hold(app.UIA3d,'on'); 
    grid (app.UIA3d,'on');
    surf(app.UIA3d,(Radius+zi*ScaleF).*cos(-xi/Radius),(Radius+zi*ScaleF).*sin(-xi/Radius),yi,zi, 'linestyle','none');
    capindex=2
    if capindex==1 %apply cap on top
    [T,R] = meshgrid(linspace(0,2*pi,64),linspace(0,app.Radius,16));
    X = R.*cos(T);
    Y = R.*sin(T);
    Z = X.*0+Y.*0+app.Height;
    dface=X.*0+Y.*0+0;
    surf(app.UIA3d,X,Y,Z,dface,'LineStyle',"none")
    end
    xlabel(app.UIA3d,'X (in)'), ylabel(app.UIA3d,'Y (in)'), zlabel(app.UIA3d,'Height (in)');
    txt0 ={'0^{o}'};
    text(app.UIA3d,0.8*Radius,0,app.Height/2,txt0,'Color','red','FontSize',13); %,'HorizontalAlignment','left')    
    txt01 ={'90^{o}'};
    text(app.UIA3d,0,-0.8*Radius,0, txt01,'Color','red','FontSize',13);    
    txt00 ={'180^{o}'};
    text(app.UIA3d,-Radius,0,txt00,'Color','red','FontSize',13);        
    txt01 ={'270^{o}'};
    text(app.UIA3d,0, 0.8*Radius,txt01,'Color','red','FontSize',13);    
    light(app.UIA3d);
    colorbar (app.UIA3d);
    rotate3d(app.UIA3d,'on');
    hold(app.UIA3d,'off');
elseif app.heightflag==1    %2D Plot at specific measuring plane
    app.heightflag=2;
    if app.partflag   %1: already have part number; 2: not 
        cla(app.UIA2d);
        theta1=[0:2*pi/360:2*pi];
        hei=app.HeightEdit;
        hold(app.UIA2d,'on');
        grid (app.UIA2d,'on');
        plot(app.UIA2d,Radius*cos(-theta1),Radius*sin(-theta1),'LineWidth',1,'Color','r','LineStyle',':') % plot nomial
        if hei==0 % plot all measured data
            i_plane=1;
            while i_plane< app.NP+1
                plot(app.UIA2d,app.LG_Scaled(i_plane,:).*cos(-theta),app.LG_Scaled(i_plane,:).*sin(-theta),'LineWidth',1,'LineStyle','-'); % plot all measurement data 
                i_plane=i_plane+1;
            end
            legendCell =['Nomial';cellstr(num2str(round(app.planeheight,2), 'h=%.2f'))];
            txt = {['Measurement:'],['Mean=' num2str(mean2(LGD))],['Max=' num2str(max(max(LGD)))],['Min=' num2str(min(min(LGD)))],['STD=' num2str(std2(LGD))]};
        else
            iplane = find(app.HeightEdit==app.planeheight);
            LG_Scaled_h = app.LG_Scaled(iplane,:);
            plot(app.UIA2d, LG_Scaled_h.*cos(-theta),LG_Scaled_h.*sin(-theta) ,'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',5) % measured at iplane
            ih= round (hei/(app.Height/double(nh)));
            plot(app.UIA2d, (Radius+zi(ih,:)*ScaleF).*cos(-xi(ih,:)/Radius),(Radius+zi(ih,:)*ScaleF).*sin(-xi(ih,:)/Radius),'Color','b','LineStyle','-'); % fitted at iplane
            legendCell ={'Nomial','Measured','Fitted'};
            txt = {['Fitted:'],['Mean=' num2str(mean(zi(ih,:)))],['Max=' num2str(max(zi(ih,:)))],['Min=' num2str(min(zi(ih,:)))],['STD=' num2str(std(zi(ih,:)))]};
        end
        grid (app.UIA2d, 'on');
        for i=1:12          
            txt0 ={(i-1)*30};
            text(app.UIA2d,0.9*Radius*cos(-(i-1)*pi/6),0.9*Radius*sin(-(i-1)*pi/6),txt0,'Rotation',180-(i-1)*30,'Color','red','FontSize',10);            
        end
        if hei==0 % plot all measured data
            legend (app.UIA2d,legendCell,'Location','southeastoutside','color','none','fontSize',6); %'Orientation','horizontal'
        else
            legend (app.UIA2d,legendCell,'Location','eastoutside','color','none','fontSize',6); %'Orientation','horizontal'
        end
        text(app.UIA2d,-app.Radius/2.2,0.1*app.Radius,txt,'FontSize',10); %,'HorizontalAlignment','left')        
        hold(app.UIA2d,'off');
    end %end of partflag
end    %end of function - mainplot ()
   
end % end of methods

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
%             drawnow; % for speeding up, not helpful
%             app.CorningIncLaserGaugeDataAnalysisUIFigure.WindowState = 'normal'; % 'fullscreen'; % 'maximized'; %screen size, not needed for WebApp    
        end

        % Callback function
        function ButtonPushed(app, event)
            app.PartName = uidropdown(app.PartID);
        end

        % Callback function
        function filereadButtonPushed(app, event)

        end

        % Button pushed function: ImportdatafileButton
        function ImportdatafileButtonPushed2(app, event)
            Readfile_ExtractData (app,event)
        end

        % Button pushed function: plot
        function plotButtonPushed(app, event)
            cla(app.UIA2d);
            cla(app.UIA3d);
            cla(app.UIA2dcontour);
            if app.dataflag==2
                display ('You must import data file first');
                return;
            else
                mainplot (app,event)
            end
        end

        % Callback function
        function PartNameValueChanged(app, event)
                    
        end

        % Value changed function: Scalefactor
        function ScalefactorValueChanged(app, event)
            value = app.Scalefactor.Value;            
        end

        % Value changed function: PartHeight
        function PartRadiusValueChanged(app, event)
        
        end

        % Value changed function: PartRadius
        function PartRadiusValueChanged2(app, event)
                                
        end

        % Value changed function: PlaneheightDropDown
        function HeightinEditValueChanged(app, event)
            if app.partflag==2
                display ('You must choose a Part ID and then run "Plot" first');
                return;
            else
                if app.PlaneheightDropDown.Value == "All"
                    app.HeightEdit =0;
                else
                    app.HeightEdit = app.planeheight(str2num(app.PlaneheightDropDown.Value(2:2))); % str2num(app.PlaneheightDropDown.Value(2:2));
                end
                app.heightflag=1;
                mainplot(app,event);
            end
        end

        % Close request function: CorningIncLaserGaugeDataAnalysisUIFigure
        function CorningIncLaserGaugeDataAnalysisUIFigureCloseRequest(app, event)
            delete(app)            
        end

        % Button pushed function: ExportresultsButton
        function ExportresultsButtonPushed(app, event)
        app.Label2.Text = "Please waite while processing your data..."
        pause(0.5)
        resultsexport(app,event);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create CorningIncLaserGaugeDataAnalysisUIFigure and hide until all components are created
            app.CorningIncLaserGaugeDataAnalysisUIFigure = uifigure('Visible', 'off');
            app.CorningIncLaserGaugeDataAnalysisUIFigure.Position = [100 100 649 491];
            app.CorningIncLaserGaugeDataAnalysisUIFigure.Name = 'Corning Inc. Laser Gauge Data Analysis';
            app.CorningIncLaserGaugeDataAnalysisUIFigure.CloseRequestFcn = createCallbackFcn(app, @CorningIncLaserGaugeDataAnalysisUIFigureCloseRequest, true);

            % Create UIA2dcontour
            app.UIA2dcontour = uiaxes(app.CorningIncLaserGaugeDataAnalysisUIFigure);
            xlabel(app.UIA2dcontour, 'Unrolled clockwise circumference in degree ({^o})')
            ylabel(app.UIA2dcontour, 'Height (in)')
            app.UIA2dcontour.XTickLabelRotation = 0;
            app.UIA2dcontour.YTickLabelRotation = 0;
            app.UIA2dcontour.ZTickLabelRotation = 0;
            app.UIA2dcontour.Position = [18 8 614 163];

            % Create UIA2d
            app.UIA2d = uiaxes(app.CorningIncLaserGaugeDataAnalysisUIFigure);
            title(app.UIA2d, 'Radius deviation - 2D')
            xlabel(app.UIA2d, 'X (in)')
            ylabel(app.UIA2d, 'Y (in)')
            app.UIA2d.XTickLabelRotation = 0;
            app.UIA2d.YTickLabelRotation = 0;
            app.UIA2d.ZTickLabelRotation = 0;
            app.UIA2d.FontSize = 11;
            app.UIA2d.Position = [323 174 327 256];

            % Create UIA3d
            app.UIA3d = uiaxes(app.CorningIncLaserGaugeDataAnalysisUIFigure);
            title(app.UIA3d, 'Surface deviation-3D')
            xlabel(app.UIA3d, 'X-Y plane')
            ylabel(app.UIA3d, 'Axial height')
            app.UIA3d.XTickLabelRotation = 0;
            app.UIA3d.YTickLabelRotation = 0;
            app.UIA3d.ZTickLabelRotation = 0;
            app.UIA3d.FontSize = 11;
            app.UIA3d.Position = [8 174 309 256];

            % Create RadiusEditFieldLabel
            app.RadiusEditFieldLabel = uilabel(app.CorningIncLaserGaugeDataAnalysisUIFigure);
            app.RadiusEditFieldLabel.HorizontalAlignment = 'right';
            app.RadiusEditFieldLabel.Position = [3 462 43 22];
            app.RadiusEditFieldLabel.Text = 'Radius';

            % Create PartRadius
            app.PartRadius = uieditfield(app.CorningIncLaserGaugeDataAnalysisUIFigure, 'numeric');
            app.PartRadius.Limits = [0.1 Inf];
            app.PartRadius.ValueChangedFcn = createCallbackFcn(app, @PartRadiusValueChanged2, true);
            app.PartRadius.Position = [49 462 36 22];
            app.PartRadius.Value = 3.22;

            % Create HeightEditFieldLabel
            app.HeightEditFieldLabel = uilabel(app.CorningIncLaserGaugeDataAnalysisUIFigure);
            app.HeightEditFieldLabel.HorizontalAlignment = 'right';
            app.HeightEditFieldLabel.Position = [84 462 40 22];
            app.HeightEditFieldLabel.Text = 'Height';

            % Create PartHeight
            app.PartHeight = uieditfield(app.CorningIncLaserGaugeDataAnalysisUIFigure, 'numeric');
            app.PartHeight.Limits = [0.1 Inf];
            app.PartHeight.ValueChangedFcn = createCallbackFcn(app, @PartRadiusValueChanged, true);
            app.PartHeight.Position = [128 462 28 22];
            app.PartHeight.Value = 5.5;

            % Create ScalefactorEditFieldLabel
            app.ScalefactorEditFieldLabel = uilabel(app.CorningIncLaserGaugeDataAnalysisUIFigure);
            app.ScalefactorEditFieldLabel.HorizontalAlignment = 'right';
            app.ScalefactorEditFieldLabel.Position = [156 462 69 22];
            app.ScalefactorEditFieldLabel.Text = 'Scale factor';

            % Create Scalefactor
            app.Scalefactor = uieditfield(app.CorningIncLaserGaugeDataAnalysisUIFigure, 'numeric');
            app.Scalefactor.Limits = [0 100];
            app.Scalefactor.ValueChangedFcn = createCallbackFcn(app, @ScalefactorValueChanged, true);
            app.Scalefactor.Position = [228 462 26 22];
            app.Scalefactor.Value = 3;

            % Create ImportdatafileButton
            app.ImportdatafileButton = uibutton(app.CorningIncLaserGaugeDataAnalysisUIFigure, 'push');
            app.ImportdatafileButton.ButtonPushedFcn = createCallbackFcn(app, @ImportdatafileButtonPushed2, true);
            app.ImportdatafileButton.BackgroundColor = [0 1 1];
            app.ImportdatafileButton.Position = [341 462 100 22];
            app.ImportdatafileButton.Text = 'Import data file';

            % Create plot
            app.plot = uibutton(app.CorningIncLaserGaugeDataAnalysisUIFigure, 'push');
            app.plot.ButtonPushedFcn = createCallbackFcn(app, @plotButtonPushed, true);
            app.plot.BackgroundColor = [0 1 1];
            app.plot.Position = [154 433 40 22];
            app.plot.Text = 'Plot';

            % Create PartIDLabel
            app.PartIDLabel = uilabel(app.CorningIncLaserGaugeDataAnalysisUIFigure);
            app.PartIDLabel.HorizontalAlignment = 'right';
            app.PartIDLabel.Position = [18 433 43 22];
            app.PartIDLabel.Text = 'Part ID';

            % Create PartName
            app.PartName = uidropdown(app.CorningIncLaserGaugeDataAnalysisUIFigure);
            app.PartName.Items = {'Part ID'};
            app.PartName.Position = [76 433 71 22];
            app.PartName.Value = 'Part ID';

            % Create PlaneheightDropDownLabel
            app.PlaneheightDropDownLabel = uilabel(app.CorningIncLaserGaugeDataAnalysisUIFigure);
            app.PlaneheightDropDownLabel.HorizontalAlignment = 'right';
            app.PlaneheightDropDownLabel.Position = [569 408 72 22];
            app.PlaneheightDropDownLabel.Text = 'Plane height';

            % Create PlaneheightDropDown
            app.PlaneheightDropDown = uidropdown(app.CorningIncLaserGaugeDataAnalysisUIFigure);
            app.PlaneheightDropDown.ValueChangedFcn = createCallbackFcn(app, @HeightinEditValueChanged, true);
            app.PlaneheightDropDown.Position = [582 386 59 22];

            % Create Label
            app.Label = uilabel(app.CorningIncLaserGaugeDataAnalysisUIFigure);
            app.Label.Position = [452 462 189 22];
            app.Label.Text = '';

            % Create ExportresultsButton
            app.ExportresultsButton = uibutton(app.CorningIncLaserGaugeDataAnalysisUIFigure, 'push');
            app.ExportresultsButton.ButtonPushedFcn = createCallbackFcn(app, @ExportresultsButtonPushed, true);
            app.ExportresultsButton.BackgroundColor = [0.8 0.8 0.8];
            app.ExportresultsButton.Position = [212 433 92 22];
            app.ExportresultsButton.Text = 'Export results';

            % Create Label2
            app.Label2 = uilabel(app.CorningIncLaserGaugeDataAnalysisUIFigure);
            app.Label2.FontSize = 14;
            app.Label2.FontWeight = 'bold';
            app.Label2.FontAngle = 'italic';
            app.Label2.FontColor = [1 0 0];
            app.Label2.Position = [307 433 281 22];
            app.Label2.Text = '';

            % Create AngleLabel
            app.AngleLabel = uilabel(app.CorningIncLaserGaugeDataAnalysisUIFigure);
            app.AngleLabel.HorizontalAlignment = 'right';
            app.AngleLabel.Position = [260 462 46 22];
            app.AngleLabel.Text = 'Angle #';

            % Create AngleNo
            app.AngleNo = uieditfield(app.CorningIncLaserGaugeDataAnalysisUIFigure, 'numeric');
            app.AngleNo.Limits = [0.1 Inf];
            app.AngleNo.Position = [309 462 25 22];
            app.AngleNo.Value = 24;

            % Show the figure after all components are created
            app.CorningIncLaserGaugeDataAnalysisUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Laser_gauge_data_analysis_v3_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.CorningIncLaserGaugeDataAnalysisUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.CorningIncLaserGaugeDataAnalysisUIFigure)
        end
    end
end