%% CP_Scan_McGrathPerrin2022.m
% This code is meant to perform line scans on the cuticular plates of
% auditory hair cells selected by the user. 
%
% (1) The ROIs are selected by clicking the center of hair cells, region 
% is duplicated, permuted, and interpolated to mimick an ImageJ reslice.
% HINT: Note the wI and hI parameters (search region width and height) when
% making selections. Try to get the bundle and a clear region in front of
% the bundle in your selection.
%
% (2) The ROIs are then analyzed. (a) The thickness and intensity of the CP
% relative to stereocilia are determined. (b) Two edges of the CP are
% determined by taking the thickness cutoffs and fitting them with a linear
% polynomial equation. The start and stop of the line scans are taken as
% two other edges. This forms a parallelogram for which height can be
% calculated and used to approximate CP thickness while accounting for the
% imaging angle.

function im = CP_Scan_McGrath2022(image_path, CellType, answer, XY_res, Z_res, UW, UH)

%% (0) Example inputs
% image_path = 'G:\Folder\File.tif';
% CellType = 'OHC';
% answer{1} = 'Animal identifier';
% answer{2} = 'Animal condition'; % Such as, 'C' for control or 'M' for mutant
% XY_res = 0.0426727; % um/px
% Z_res = 0.1502039; % um/px
% UW = 20; % 20 pixels * XY_res in nm/px ~ 0.4um
% UH = 100; % 250 pixels * XY_res in nm/px ~ 10um

%% (1) Select ROIs in image and save data to .mat output.
% Get image info and user selections
[path, name, ext] = fileparts(image_path);
file = [name{:} ext{:}];

% Import the tiff stacks and place info into "im" structure
im.Path = fullfile(path, file);
Data3D = FastTiff(im.Path);
im.Data2D = max(Data3D,[],3);
Info = imfinfo(im.Path, 'tif');
[im.Width, w] = deal(Info(1).Width); 
[im.Height, h] = deal(Info(1).Height);

% Call figure and select the centers of the hair cell apical surface. This
% will be used to crop out a "w" by "h" square section out of the stack,
% including all Z-slices in the crop. These images will be used to generate
% the figures shown to the user to select CP and stereocilia rois.
f = figure; f.Position = [100 300 w h];
imagesc(im.Data2D); colormap('gray');
title(['Select the ' CellType])
[im.xi, im.yi] = getpts(f);
close(f);

% Crop out selections
I = length(im.xi); Crop3D = cell(I, 1);
for i = 1:I
    x1 = ceil(im.xi(i) - 0.5*UW);
    x2 = ceil(im.xi(i) + 0.5*UW);
    y1 = ceil(im.yi(i) - 0.5*UH);
    y2 = ceil(im.yi(i) + 0.5*UH);
    var = [x1, x2, y1, y2];
    var(var < 0) = 0;
    Crop3D{i} = Data3D(y1:y2, x1:x2, :);
end
im.Crop3D = Crop3D;

% Reslice selections and interpolate. This is similar to ImageJ's ReSlice
% function. I used the default "spline" interpolation settings in Matlab.
clear tmp tmp2

% Save path for reslice image data
savepath = fullfile(path, "Reslices");
if ~exist(savepath,'dir')
    mkdir(savepath);
end

tic;
I = length(im.Crop3D); RSPaths = cell(I, 1);
for i = 1:I
    
    pmIm = rescale(permute(im.Crop3D{i}, [3, 1, 2]), 0, 1);
    [m, n, o] = size(pmIm);
    tmppmIm = zeros(m, n, o);
    
    sc = 2;
    xq  = 1:(XY_res/Z_res):m;
    method = 'spline';
    
    tmp = interp1(tmppmIm, xq, method);
    m2 = 2^sc * (size(tmp, 1) - 1) + 1;
    n2 = 2^sc * (size(tmp, 2) - 1) + 1;
    tmp2 = zeros(m2, n2, o);
    
    % First interpolate to account for poor Z-resolution (square voxels to
    % XY_res nm^3)
    for j = 1:n
        for k = 1:o
            v = pmIm(:, j, k); 
            tmp(:, j, k) = interp1(v, xq, method);
        end
    end
    
    % Next interpolate each frame to smooth the image appearance
    for j = 1:o
        v = tmp(:, :, j);
        tmp2(:, :, j) = interp2(v, sc, method);
    end
    name = ['Cell_' num2str(i) '.mat'];
    imdata = tmp2;
    RSPaths{i} = fullfile(savepath, name);
    save(RSPaths{i}, 'imdata');
    
end
disp(['Reslice complete and saved (' num2str(toc) 's)'])

im.RSPaths = RSPaths;

%% (2) Get user selection of stereocilia (rectangle) and cuticular plate
% (line scans) ROIs

locbuff = 0.5; 
% ^^^ This is the fraction of the most common pixel value in the
% cuticular plate that will be used to define it's edges. If it is set to 
% 0.7 then the pixel at 70% the common value will define the edge. This
% tends to be below the peaks and troughs caused by pixel variance in the 
% phalloidin labeling signal. Less variance allows for tighter cutoffs,
% more variance allows for looser cutoffs.

% Select the stereocilia ROI
I = length(RSPaths); 
[~, SCint, ~] = deal(cell(I, 1));
for i = 1:I
    
    matObj = matfile(RSPaths{i});
    RS = matObj.imdata;
    
    % Make figure showing the reslice and await instruction to continue.
    fig = figure;
    imshow(rescale(sum(RS, 3)))
    colormap('turbo')
    f = msgbox('Click "OK" after selecting ROI.');
    title('Select the stereocilia bundle.','FontSize',18)
    set(fig,'Position',[91   173   841   574]);
    set(f,'Position',[741.0000  320.2500  150.0000   51.7500]);
    A = drawrectangle('Color','m');
    uiwait(f); 
    
    % Crop out the selection. Surely there is a better way. I didn't use
    % it.
    im.B0{i} = A.Position;
    B = im.B0{i};
    width = B(3);
    height = B(4);
    c = images.spatialref.Cuboid(...
        [B(1),B(1)+width],[B(2),B(2)+height],[1,size(RS,3)]);
    ROI =...
        imcrop3(RS,c);
    logicA = rescale(ROI, 0, 1) < 0.90;
    tmpROI = ROI;
    tmpROI(logicA) = nan; 
    SCint{i} = tmpROI; %clear logicA tmpROI
    
    disp(['Stereocilia ' num2str(i) '/' num2str(I) ' ROI recorded']);
    close(fig);
    im.ROI{i} = ROI;
    
end

im.SCint = SCint;

% Select the cuticular plate ROI
for i = 1:I
    
    matObj = matfile(RSPaths{i});
    RS = matObj.imdata;
    
    fig = figure;
    imshow(rescale(sum(RS, 3)))
    colormap('turbo')
    f = msgbox('Click "OK" after selecting ROI.');
    title('Draw line perpendicularly across the cuticular plate.','FontSize',18)
    set(fig,'Position',[91   173   841   574]);
    set(f,'Position',[741.0000  320.2500  150.0000   51.7500]);
    ROI = drawline;
    uiwait(f); 
    LINE{i} = ROI.Position; 
    disp('Cuticular plate ROI recorded');
    close(fig);
    
end

im.LINE = LINE;

%% (3) Show and compile results
% Get image scaling and line lengths
L = table(...
    zeros(I, 1), zeros(I, 1), zeros(I, 1), zeros(I, 1), zeros(I, 1),...
    'VariableNames', {'x1', 'x2', 'y1', 'y2', 'Length'});
for i = 1:I
    
    matObj = matfile(RSPaths{i});
    RS = matObj.imdata;
    
    % Get scaling factors
    [M1, ~, O1] = size(im.Crop3D{1});
    [O2, M2, ~] = size(RS);
    ySF = M1/M2;
    xSF = O1/O2;
    
    A = LINE{i}; % Get cuticular plate ROI
    L.x1(i) = A(1); L.x2(i) = A(2);
    L.y1(i) = A(3); L.y2(i) = A(4);
    
    % Determine line lengths
    A = L(i, :);
    L.Length(i) = sqrt(((A.x1 - A.x2)*xSF)^2 + ((A.y1 - A.y2)*ySF)^2);
    
end

% Process cuticular plate and junction data and show user readout
for i = 1:I
    
    matObj = matfile(RSPaths{i});
    RS = matObj.imdata;
    
    fig = figure(i);
    T = tiledlayout(1, 2);
    TITLE = ['Organism ' answer{1}...
    ' - Condtion ' answer{2} ' - Cell ' num2str(i)];
    title(T, TITLE)
    
    nexttile
    A = L(i,:);
    
    title('Line scan on image')
    imshow(rescale(sum(RS, 3)), 'InitialMagnification', 400)
    colormap('turbo')
    hold on
    plot([A.x1, A.x2], [A.y1, A.y2], 'c-', 'LineWidth', 3)
    hold off
    
    % Process data, plot "cuticular plate" profile, and annotate
    nexttile
    title('Cuticular plate profile')
    ylabel('Norm. Intensity')
    xlabel('Distance (nm)')
    A = L(i,:);
    
    hold on
    [~, ~, J] = size(RS); 
    map = turbo(J); 
    [ypos1, ypos2, LineScan] = deal(cell(J, 1));
    [thickness, CPintTmp] = deal(zeros(J, 1));
    clear tmpLineScan
    for j = 1:J
        c = improfile(RS(:, :, j),...
        [A.x1 A.x2], [A.y1 A.y2]);
        clear LineScan
        LineScan(:, 1) = (rescale(1:length(c)).*A.Length).*XY_res;
        LineScan(:, 2) = rescale(c);
        LineScan(:, 3) = c;
        [f, xi] = ksdensity(interp1(LineScan(:, 2),1:0.5:length(LineScan(:, 2))));
        [pks, loc] = findpeaks(f, xi, 'SortStr', 'descend');
        [~, idx] = max(loc);
        
        % Take the pixels from the 3D selection that pass the criteria and
        % get mean intensity from those pixels
        tmpInt = LineScan(LineScan(:, 2) > locbuff*pks(idx), 3);
        
        % Get the CP thickness from edge to edge across the x-axis
        tmpPos = find(LineScan(:, 2) > locbuff*max(loc));
        ypos1{j} = LineScan(tmpPos(1), 1);
        ypos2{j} = LineScan(tmpPos(end), 1);
        thickness(j) = abs(ypos2{j} - ypos1{j});
        CPintTmp(j) = mean(tmpInt)/mean(SCint{i}, 'all', 'omitnan');
        
        plot(LineScan(:, 1), LineScan(:, 3), 'Color', map(j, :))
        xline(ypos1{j}, '--', 'Color', map(j, :))
        xline(ypos2{j}, '--', 'Color', map(j, :))
        tmpLineScan{j} = LineScan;
    end
    hold off
    im.CP.Thickness(i) = mean(thickness);
    im.CP.ThicknessSTD(i)  = std(thickness);
    im.CP.LineScan{i} = tmpLineScan;
    im.CP.Intensity(i) = mean(CPintTmp, 'all', 'omitnan');
    
    CellPath = fullfile(path, 'Data pics');
    if ~exist(CellPath, 'dir')
        mkdir(CellPath)
    end
    bdataFolder = fullfile(path, 'Bulk data');
    if ~exist(bdataFolder, 'dir')
        mkdir(bdataFolder);
    end
    
    CellName = [answer{1} '_' CellType '_' num2str(i) '.pdf'];
    dataName = [answer{1} '_' CellType '_' num2str(i)];
    save(fullfile(bdataFolder, [dataName '_' 'thickness' '.mat']), 'thickness');
    save(fullfile(bdataFolder, [dataName '_' 'intensity' '.mat']), 'CPintTmp'); 
    save(fullfile(bdataFolder, [dataName '_' 'linescan' '.mat']), 'tmpLineScan'); 
    
    Grab = gcf;
    set(Grab, 'Units', 'normalized');
    set(Grab, 'Position', [0 0 1 1]);
    set(Grab, 'PaperUnits', 'normalized');
    set(Grab, 'PaperPosition', [0 0 1 1]);
    set(Grab, 'PaperOrientation', 'landscape');
    print(fig, fullfile(CellPath,CellName), '-dpdf')
    close(fig)
    
end

save(fullfile(bdataFolder,[answer{1} CellType '_data.mat']), 'im');

%% Function(s)

function data = FastTiff(filename)
    warning('off','all') % Suppress all the tiff warnings
    tstack  = Tiff(filename);
    [II,JJ] = size(tstack.read());
    KK = length(imfinfo(filename));
    data = zeros(II,JJ,KK);
    data(:,:,1)  = tstack.read();
    for nn = 2:KK
        tstack.nextDirectory()
        data(:,:,nn) = tstack.read();
    end
    warning('on','all')
end

%% EOF

end