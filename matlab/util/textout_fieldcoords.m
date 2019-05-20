function textout_fieldcoords()

  paths = get_paths();

  load(sprintf('%scataloginfo.mat',paths.catdir)); 
  [~,n] = size(field_number);

  [gout] = coco([field_RA',field_DEC'],'J2000.0','g','d','d');
  
  glon = gout(:,1);
  glat = gout(:,2);
  
  [eout] = coco([field_RA',field_DEC'],'J2000.0','e','d','d');
  
  elon = eout(:,1);
  elat = eout(:,2);

  fid = fopen('lookup/nh_fieldcoords.txt','w');
  for ifield=1:n
    mystring = sprintf('%d, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f\n',...
	ifield,field_RA(ifield),field_DEC(ifield),glon(ifield),glat(ifield),...
	elon(ifield),elat(ifield));
    fprintf(fid,mystring);
  end
  fclose(fid);
  
  
%end
