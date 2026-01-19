function parcel_to_HCPsurface(parcel_value,output_name,color_range,cmap,table_path)




[NUM,~,RAW] = xlsread(table_path)
Nemo_indx = NUM(:,1)
HCP_index = NUM(:,3)
parcel_value_HCP = zeros(1,360);
parcel_value_HCP(HCP_index) = parcel_value(Nemo_indx)
parcel_value_fsa5 = parcel_to_surface(parcel_value_HCP, 'glasser_360_fsa5');
parcel_value_fsa5(find( parcel_value_fsa5 == 0)) = -1000;
% Project the results on the surface brain
f = figure,
plot_cortical(parcel_value_fsa5, 'surface_name', 'fsa5', 'color_range', ...
    color_range, 'cmap', cmap)

print(gcf,'-dtiff','-r300',output_name);
%close(gcf)

end