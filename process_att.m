%{
Testing matlab isosurface() and other functions for rendering snowflake .att files. 
utilty functions: readAtt.m  stl_write.m
%}
function process_att(threshold_val)
    if nargin < 1
        threshold_val = 0.995;
    end

    global threshold face_color face_alpha;
    global line_color line_width;

    %function process_att(opts)
    %    arguments
    %        opts.threshold (1,1) {mustBeNumeric} = 0.999
    %    end
    %util_dir = genpath("../BIN");  addpath(util_dir);

    % coplanarity threshold (two surface normals inner product)
    threshold = threshold_val;

    % Morphological operation
    % opt_wgt = 'imopen';  
    opt_wgt = 'imclose';  
    % 

    cran = "#00FFFF";
    lightcran = "#E0FFFF";
    azure = "#F0FFFF";  %Azure
    ghostwhite = "#F8F8FF";

    face_color = cran;
    face_alpha = 0.4;

    line_color = lightcran;  
    line_width = 0.5;


    %distances to the center
    x1 = 0:2;
    [x2,y2,z2] = meshgrid(x1,x1,x1);
    distc = (x2-x2(2,2,2)).^2+(y2-y2(2,2,2)).^2+(z2-z2(2,2,2)).^2;
    distc = distc.^0.5;
    b_equal = 1.0 ./distc;

    %Gaussian weights
    mu=0;
    sig=1.0;
    b_wgt = normpdf(distc,mu,sig);
    sumc = sum(b_wgt,'all');
    b_conv = b_wgt/sumc;

    %unit_test 
    %[so,sf_i] = unit_test(Name='cube');
    %so = unit_test(Name='o3d_bunny');
    %return

    fid = fopen('list_att_files.txt','r');
    while ~feof(fid)
        tline = fgetl(fid);
        if ismember(tline(1:1),{'#','!','%'}) 
           continue;
        end 
        dat = readAtt(tline);
        [folder, baseNoExt, extension] = fileparts(tline);

        dim = size(dat);
        if contains(opt_wgt, 'Equal', 'IgnoreCase', true)
            %averaging by 3x3x3 cube equal weighting
            dat = convn(dat,b_equal,'same');
        elseif contains(opt_wgt, 'Gauss', 'IgnoreCase', true)
            %averaging by 3x3x3 cube equal weighting
            dat = convn(dat,b_conv,'same');
        elseif contains(opt_wgt, 'Dilate', 'IgnoreCase', true)
            se_opt = 'cube';
            se = strel(se_opt, 3); 
            dat = imdilate(dat, se);
            opt_wgt = ['_Dilate' se_opt];
        elseif contains(opt_wgt, 'imopen', 'IgnoreCase', true)
            se_opt = 'cube';
    	    se = strel(se_opt, 3); 
            dat = imopen(dat, se);
            dat = convn(dat,b_conv,'same');
            opt_wgt = ['_imopen' se_opt 'Gauss'];
        elseif contains(opt_wgt, 'imclose', 'IgnoreCase', true)
            se_opt = 'cube';
    	    se = strel(se_opt, 3); 
            dat = imclose(dat, se);
            dat = convn(dat,b_conv,'same');
            opt_wgt = ['_imclose' se_opt 'Gauss'];
        elseif contains(opt_wgt, 'bwperim', 'IgnoreCase', true)
            se_opt = '26';
            dat = bwperim(dat, 26);
            opt_wgt = ['_bwperim' se_opt];
        else
            %averaging by 3x3x3 cube equal weighting
            dat = convn(dat,b_conv,'same');
        end


        dat_1 = zeros(dim+2);
        dat_1(2:dim(1)+1,2:dim(2)+1,2:dim(3)+1) = dat;

        %isosurface plot 
        iso = plot_iso_sur(dat_1,baseNoExt,opt_wgt);

        %save to stl file
        IsOutSTL = 0;
        if IsOutSTL
            fstl=[baseNoExt '.stl']; 
            stl_write(fstl, iso);

            %%Read the stl file to visualize the output
            f2 = figure(2);
            TO = stlread(fstl);
            trimesh(TO);
            daspect([1 1 1]);
            ca = [0 -1 0];
            %[az,el] = deal(0.6438,-49.3887); view(az,el);
            view(ca);
            print('-dpng','-r200','-loose', [baseNoExt '_stl'] );
            %close all;
        end
    end
end


%
function so = plot_iso_sur(dat_v,fn_key,opt_w)
    global threshold face_color face_alpha
    global line_color line_width;

    matfile = sprintf('r03_so_%s_%s%s.mat', fn_key,num2str(threshold),opt_w);
    matfile_toLoad = sprintf('r03_so_%s_%s%s_toLoad.mat', fn_key,num2str(threshold),opt_w);
    if exist(matfile_toLoad, 'file')
        ma = load(matfile);
        vertices = ma.vertices;
        faces = ma.faces;
        coplanarGroups = ma.coplanarGroups;
	outerEdges = ma.outerEdges;
        disp(['Coplanar face groups loaded: ', matfile_toLoad]);
    else
        [nx, ny, nz] = size(dat_v);
        [X, Y, Z] = meshgrid(1:ny, 1:nx, 1:nz);

        % Add Gaussian noise perturbation (2% of voxel size 1), clipped at 3*sigma
        sigma = 0.02 / 3;
        limit = 3 * sigma;

        X = X + max(min(sigma * randn(size(X)), limit), -limit);
        Y = Y + max(min(sigma * randn(size(Y)), limit), -limit);
        Z = Z + max(min(sigma * randn(size(Z)), limit), -limit);

        so = isosurface(X,Y,Z,dat_v,0.5);
        vertices = so.vertices;
        faces = so.faces;
	[coplanarGroups,outerEdges] = get_coplanar_groups(vertices,faces,fn_key);
        save(matfile,'coplanarGroups','outerEdges','vertices','faces');
        disp(['Coplanar face groups saved: ', matfile]);
    end

    % Create a custom colormap with shades of blue
    numGroups = length(coplanarGroups);
    % Create a gradient from dark blue to bright blue
    blueShades = [zeros(numGroups, 1), linspace(0.1, 0.5, numGroups)', linspace(0.3, 1, numGroups)'];

    % Initialize the color array for all faces
    numFaces = size(faces, 1);
    faceColorData = zeros(numFaces, 3);
    
    % Assign the appropriate color to each face based on its group
    for k = 1:numGroups
        groupIndices = coplanarGroups{k};
        groupColor = blueShades(k, :);
        faceColorData(groupIndices, :) = repmat(groupColor, length(groupIndices), 1);
    end

    set(0,'DefaultFigureVisible','off')
    figure; hold on; grid on;
    %daspect([1,1,1]);
    axis equal;
    %ca = [0 -1 0]; view(ca); 
    view(3);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    %camlight
    camlight(-80,-10);
    lighting gouraud
   
    % Display the patch for the current group with no internal edges
    so = patch('Faces', faces, 'Vertices', vertices, ...
          'FaceColor', face_color,   ... %'FaceVertexCData', faceColorData, ... % Assign the color map
          'FaceAlpha', face_alpha, ...
          'EdgeColor', 'none');
    %gso.FaceAlpha = 0.5;
    
    % Loop through each identified group to plot it
    for k = 1:numGroups
        boundaryEdges = outerEdges{k}; 
    
        if ~isempty(boundaryEdges)
            edgeVertices = vertices(boundaryEdges', :);
            X = reshape(edgeVertices(:,1), 2, []);
            Y = reshape(edgeVertices(:,2), 2, []);
            Z = reshape(edgeVertices(:,3), 2, []);
            %plot3(X, Y, Z, 'b-', 'LineWidth', 0.5); 
            plot3(X, Y, Z, '-','Color',line_color, 'LineWidth', line_width); 
        end
    end
    
    hold off;
    
    disp(['Identified ', num2str(numGroups), ' planar groups.']);
 
    ftag = sprintf('Particle: %s,  threshold=%s, Coplanar=%s', ...
           strrep(fn_key, '_', '\_'), num2str(threshold),num2str(numGroups));
    title(ftag, 'FontSize',18);

    cf = gcf; ax = gca;
    ax.FontSize = 14; 
    set(cf,'Outerposition',[1 1 1600 1000]);
    pngfile = sprintf('r04_so_%s_%s%s.png', fn_key,num2str(threshold),opt_w);
    print('-dpng','-r800','-loose', pngfile );
    disp(['png file generated: ', pngfile]);
    figfile = sprintf('r05_so_%s_%s%s.fig', fn_key,num2str(threshold),opt_w);
    savefig(figfile);
    set(0,'DefaultFigureVisible','on')
end
%
function [coplanarGroups,outerEdges] = get_coplanar_groups(vertices,faces,fn_key)
%   check if two connected faces on same plane by surface normal
    global threshold
    
    warning('off','MATLAB:triangulation:PtsNotInTriWarnId');

    % Create a main triangulation object for the entire model
    mainTri = triangulation(faces, vertices);
    
    % Calculate face normals and identify neighbors
    faceNormals = faceNormal(mainTri);
    neighboringFaces = neighbors(mainTri);
    
    % Group all coplanar faces
    numFaces = size(faces, 1);
    visitedFaces = false(numFaces, 1);
    coplanarGroups = {}; % Initialize a cell array to store the groups
    
    for i = 1:numFaces
        if ~visitedFaces(i)
            currentGroupIndices = i;
            visitedFaces(i) = true;
            
            % Queue for a breadth-first search
            queue = i;
            
            while ~isempty(queue)
                currentFace = queue(1);
                queue(1) = [];
                
                currentNeighbors = neighboringFaces(currentFace, :);
                
                for j = 1:length(currentNeighbors)
                    neighbor = currentNeighbors(j);
                    
                    if ~isnan(neighbor) && ~visitedFaces(neighbor)
                        dotProduct = dot(faceNormals(currentFace,:), faceNormals(neighbor,:));
                        
                        if abs(abs(dotProduct)) > threshold % Coplanarity tolerance
                            currentGroupIndices = [currentGroupIndices; neighbor];
                            visitedFaces(neighbor) = true;
                            queue = [queue; neighbor]; 
                        end
                    end
                end
            end
            
            % Add the completed group to cell array
            coplanarGroups{end+1} = currentGroupIndices;
        end
    end

    % Outer edge for each group
    numGroups = length(coplanarGroups);
    outerEdges = {};
    for k = 1:numGroups
        currentGroupIndices = coplanarGroups{k};
        groupFaces = faces(currentGroupIndices, :);
    
        % Create a temporary triangulation for the group to find its boundary
        groupTri = triangulation(groupFaces, vertices);
	outerEdges{k} = freeBoundary(groupTri);
    end
end 

function so = unit_test(opts)
% test functions with a cube
arguments
    opts.Name (1,:) char {mustBeMember(opts.Name,{'cube','o3d_bunny'})} = 'o3d_bunny'
end
    global threshold; 
    display(opts.Name)
    if strcmp(opts.Name, 'cube')
        vert = [0 0 0;1 0 0;1 1 0;0 1 0;0 0 1;1 0 1;1 1 1;0 1 1];
        fac = [1 2 6; 6 5 1;2 3 7;7 6 2;3 4 8;8 7 3;4 1 5;5 8 4;1 2 3;3 4 1;8 7 6;6 5 8];
        fpng = 'cube_patch.png'
    elseif opts.Name(1,:) == 'o3d_bunny'
        %fn = "bunny_ms_alpha_0.135721.stl"
        %fn = "bunny_ms_alpha_0.500000.stl"
        fn = "bunny_ms_alpha_0.036840.stl"
        %fn = "bunny_ms_alpha_0.010000.stl"
        [pathstr,name,ext] = fileparts(fn);
        %%Read the stl file to visualize the output
        bms =  stlread(fn);

        %f2 = figure(2);
        %trimesh(bms);
        %daspect([1 1 1]);
        %ca = [0 -1 0]; view(ca); %view(3);
        %print('-dpng','-r200','-loose', name+'_stl.png' );
        %close all;

        vert = bms.Points;
        fac = bms.ConnectivityList;
        fpng = name+'_patch.png'
    end
    nf = size(fac,1);
    colormap(jet(128));
    colors=zeros(nf,3);
    colors(:,3) = 0.8;

    s =  struct('Vertices',vert,'Faces',fac);
    so = check_on_same_subset(s);

    so = patch(so,'FaceColor','flat',...
        'EdgeColor','#80B3FF','EdgeAlpha',0.4);
    
    daspect([1,1,1]);
    ca = [0 -1 0]; view(ca); view(3);
    grid on; set(gcf,'Outerposition',[1 1 1600 1000]);
    print('-dpng','-r200','-loose', fpng );

end
