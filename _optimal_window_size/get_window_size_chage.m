fpIn         = '/main/calen/occluding/images/filtered_images/';
bin         = 1;
eccentricity = 1:5;
targets      = {'vertical', 'horizontal', 'bowtie', 'spot'};
f            = @(x) occluding_model.lib.getTemplateResponse(Settings,fpIn,bin,eccentricity,targets,1:500,x,psf4mm)

winVals = 0:.5:2;
dat = arrayfun(f, winVals);

g = @(x) occluding_model.lib.writeTemplateMat2R(dat(x), 1:500, numel(bin), numel(targets), numel(eccentricity), ['~/Dropbox/Calen/Dropbox/', replace(num2str(winVals(x)), '.', ''), '.txt'])

arrayfun(g, 1:length(winVals));


