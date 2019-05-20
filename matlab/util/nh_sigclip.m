function thismean=nh_sigclip(indat);

  for iter=1:3
    sigma = std(indat);
    if iter==1
      mu = median(indat);
      newdat = indat;
    else
      mu = mean(newdat);
    end
    whpl = abs(newdat) < mu + 3.*sigma; 
    newdat = newdat(whpl);
  end

  thismean = mean(newdat);


% end