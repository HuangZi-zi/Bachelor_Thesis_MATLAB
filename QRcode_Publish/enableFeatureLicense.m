function enableFeatureLicense(feature)
  % by Brad Hieb
  license('test', feature, 'enable'); 
  license('checkout', feature, 'enable');
end

