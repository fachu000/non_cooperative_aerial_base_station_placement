function h = plotMyGraph(weightMatrix, tags, pipes_color, b_onlyCircles)

% Si se barajan (no barajean*) las tags, deber?an seguir saliendo en
% las mismas coordenadas xy.
N = size(weightMatrix,1);
file_name_prefix = 'vannsystem_dagens_';
file_extension = '.png';
selected14Tags = string({'AAS.FT1.SISTIM', 'BRE.FT1.SISTIM', 'BRV.FT1.SISTIM', 'GAN.FT1.SISTIM', 'HEL.FT1.SISTIM', 'KOK.FT1.SISTIM', 'KVE.FT1.SISTIM', 'OVR.FT1.SISTIM', 'SLE.FT1.SISTIM', 'SOD.FT1.SISTIM', 'TOR.FT1.SISTIM', 'EID.FT1.SISTIM','FLE.FT1.SISTIM','KRO.FT1.SISTIM'});
assert(N<=14, 'too many tags!');

% AAS.FT1.SISTIM
xCoordinate(1) = 570; yCoordinate(1) = 980;
% BRE.FT1.SISTIM
xCoordinate(2) = 470; yCoordinate(2) = 800;
% BRV.FT1.SISTIM
xCoordinate(3) = 600; yCoordinate(3) = 1100;
% GAN.FT1.SISTIM
xCoordinate(4) = 600; yCoordinate(4) = 400;
% HEL.FT1.SISTIM
xCoordinate(5) = 495; yCoordinate(5) = 825;
% KOK.FT1.SISTIM
xCoordinate(6) = 800; yCoordinate(6) = 550;
% KVE.FT1.SISTIM
xCoordinate(7) = 1050; yCoordinate(7) = 300;
% OVR.FT1.SISTIM
xCoordinate(8) = 640; yCoordinate(8) = 950;
% SLE.FT1.SISTIM
xCoordinate(9) = 600; yCoordinate(9) = 900;
% SOD.FT1.SISTIM
xCoordinate(10) = 775; yCoordinate(10) = 500;
% TOR.FT1.SISTIM
xCoordinate(11) = 1250; yCoordinate(11) = 725;
% EID.FT1.SISTIM
xCoordinate(12) = 830; yCoordinate(12) = 460;
% FLE.FT1.SISTIM
xCoordinate(13) = 900; yCoordinate(13) = 1400;
% KRO.FT1.SISTIM
xCoordinate(14) = 600; yCoordinate(14) = 800;

%x = zeros(N,1); y = zeros(N,1);

for i=1:N
    for j=1:14
        if(tags(i) == selected14Tags(j))
            x(i) = xCoordinate(j);
            y(i) = yCoordinate(j);
            break;
        end
    end
end
image_name_complete = char(strcat(file_name_prefix,pipes_color,file_extension));
img=imread(image_name_complete);

if b_onlyCircles
	
	imagesc(img);
hold on

G1=digraph(zeros(N),'OmitSelfLoops');
%h=plot(G1,'XData',xCoordinate,'YData',yCoordinate,'ArrowSize',15,'EdgeCData',G1.Edges.Weight,'LineWidth',4, 'NodeLabel', cellstr(tags));

	h=plot(G1,'XData',x,'YData',y,'ArrowSize',15, ...
		'EdgeCData',G1.Edges.Weight,'LineWidth',4, ...
		'NodeLabel', cellstr(tags), 'MarkerSize', 15, 'Marker', 'o'...
		, 'NodeColor', [1, 0.1, 0.8]);
	% gg  = gplot(weightMatrix, [x', y'], 'o');
	% set(gg, 'MarkerSize', 10);
	load cmap_arrows;
	set(gcf, 'colormap', cmap);
	
else
	
imagesc(img);
hold on

G1=digraph(weightMatrix','OmitSelfLoops');
%h=plot(G1,'XData',xCoordinate,'YData',yCoordinate,'ArrowSize',15,'EdgeCData',G1.Edges.Weight,'LineWidth',4, 'NodeLabel', cellstr(tags));

	h=plot(G1,'XData',x,'YData',y,'ArrowSize',15, ...
		'EdgeCData',G1.Edges.Weight,'LineWidth',4, ...
		'NodeLabel', cellstr(tags), 'MarkerSize', 12, 'Marker', 'o'...
		, 'NodeColor', [1, 0.1, 0.8]);
	% gg  = gplot(weightMatrix, [x', y'], 'o');
	% set(gg, 'MarkerSize', 10);
	load cmap_arrows;
	set(gcf, 'colormap', cmap);
end

