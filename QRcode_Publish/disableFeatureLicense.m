function disableFeatureLicense(feature)
  % by Brad Hieb

  % [Usage Example]      disableFeatureLicense('simulink')
  %
  % 'curve_fitting_toolbox'
  % 'distrib_computing_toolbox'   Parallel Computing Toolbox
  % 'fixed_point_toolbox'
  % 'image_acquisition_toolbox'
  % 'image_toolbox'
  % 'matlab'
  % 'signal_blocks'                DSP ST
  % 'signal_toolbox'               SIG
  % 'simulink'
  % 'statistics_toolbox'
  % 'video_and_image_blockset'     CVST
  % 'simulink_report_gen'
  % 'real-time_workshop'           Simulink Coder
  % 'rtw_embedded_coder'           Embedded Coder(旧Real-Time Workshop Embedded Coder)          Simulinkから=>'real-time_workshop'(現Simulink Coder)  'matlab_coder' ライセンスも使う
  %

  license('test', feature, 'disable');         % => license('test',feature) always return 0
  license('checkout', feature, 'disable');
end

