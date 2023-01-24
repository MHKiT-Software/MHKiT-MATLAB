classdef Dolfyn_Test_Orient < matlab.unittest.TestCase       
   
    methods (Test)
          
        function orient_testcase(testCase) 
            % These tests confirm that the euler2orient and orient2euler 
            % functions are consistent, and that they follow the 
            % conventions defined in the DOLfYN documentation 
            % (data-structure.html#heading-pitch-roll), namely:
     
            %   - a "ZYX" rotation order. That is, these variables are 
            %     computed assuming that rotation from the 
            %     earth -> instrument frame happens by rotating around the 
            %     z-axis first (heading), then rotating around the y-axis 
            %     (pitch), then rotating around the x-axis (roll).
     
            %   - heading is defined as the direction the x-axis points, 
            %     positive clockwise from North (this is the opposite 
            %     direction from the right-hand-rule around the Z-axis)
     
            %   - pitch is positive when the x-axis pitches up (this is 
            %     opposite the right-hand-rule around the Y-axis)
     
            %   - roll is positive according to the right-hand-rule 
            %     around the instument's x-axis
     
            % IF YOU MAKE CHANGES TO THESE CONVENTIONS, BE SURE TO UPDATE 
            % THE DOCUMENTATION.
            
            test_omat = [[0, 1, 0];[-1, 0, 0];[0, 0, 1]];
            test_omat = reshape(test_omat,[1,1,3,3]);
            Obj.diff = Dolfyn_Test_Orient.check_hpr(...
                0, 0, 0, test_omat);
            testCase.assertLessThan(Obj.diff, 1e-10);

            test_omat = [[1, 0, 0];[0, 1, 0];[0, 0, 1]];
            test_omat = reshape(test_omat,[1,1,3,3]);
            Obj.diff = Dolfyn_Test_Orient.check_hpr(...
                90, 0, 0, test_omat);
            testCase.assertLessThan(Obj.diff, 1e-10);

            test_omat = [[1, 0, 0];[0, 0, 1];[0, -1, 0]];
            test_omat = reshape(test_omat,[1,1,3,3]);
            Obj.diff = Dolfyn_Test_Orient.check_hpr(...
                90, 0, 90, test_omat);
            testCase.assertLessThan(Obj.diff, 1e-10);

            sq2 = 1/sqrt(2);

            test_omat = [[sq2, sq2, 0];[-sq2, sq2, 0];[0, 0, 1]];
            test_omat = reshape(test_omat,[1,1,3,3]);
            Obj.diff = Dolfyn_Test_Orient.check_hpr(...
                45, 0, 0, test_omat);
            testCase.assertLessThan(Obj.diff, 1e-10);

            test_omat = [[0, sq2, sq2];[-1, 0, 0];[0, -sq2, sq2]];
            test_omat = reshape(test_omat,[1,1,3,3]);
            Obj.diff = Dolfyn_Test_Orient.check_hpr(...
                0, 45, 0, test_omat);
            testCase.assertLessThan(Obj.diff, 1e-10);

            test_omat = [[0, 1, 0];[-sq2, 0, sq2];[sq2, 0, sq2]];
            test_omat = reshape(test_omat,[1,1,3,3]);
            Obj.diff = Dolfyn_Test_Orient.check_hpr(...
                0, 0, 45, test_omat);
            testCase.assertLessThan(Obj.diff, 1e-10);

            test_omat = [[sq2, 0, sq2];[-sq2, 0, sq2];[0, -1, 0]];
            test_omat = reshape(test_omat,[1,1,3,3]);
            Obj.diff = Dolfyn_Test_Orient.check_hpr(...
                90, 45, 90, test_omat);
            testCase.assertLessThan(Obj.diff, 1e-10);
        end

        function test_pr_declination(testCase)
            % Test to confirm that pitch and roll don't change when you set
            % declination
            declin = 15.37;

            dat = read_netcdf('../../examples/data/dolfyn/control/vector_data_imu01.nc');
            [h0, p0, r0] = orient2euler(dat);

            dat = set_declination(dat, declin);
            [h1, p1, r1] = orient2euler(dat);

            testCase.assertLessThan(abs(p1-p0), 1e-5,...
                'Pitch changes when setting declination');
            testCase.assertLessThan(abs(r1-r0), 1e-5,...
                'Roll changes when setting declination');
            testCase.assertLessThan(abs((h0+declin)-h1), 1e-5,...
                'Pitch changes when setting declination');
        end

        function test_q_hpr(testCase)
            dat = read_netcdf('../../examples/data/dolfyn/control/Sig1000_IMU.nc');
            dcm = quaternion2orient(dat.quaternions);

            diff = 0.0;
            cntrl = dat.orientmat;
            % Data
            tmp1 = isnan(cntrl.data);
            dt1 = double(cntrl.data);
            dt1(tmp1) = 0.0;
            tmp2 = isnan(dcm.data);
            dt2 = double(dcm.data);
            dt2(tmp2) = 0.0;

            diff = diff + abs(sum(abs(dt1 - dt2),...
                    1:numel(size(dt1)))/length(dt1)); 
            % Dims
            for kk = 1:length(cntrl.dims)
                diff = diff + double(~strcmpi(...
                    cntrl.dims{kk}, dcm.dims{kk}));
            end
            % Coords 
            cntrl_coord_names =fieldnames(cntrl.coords);
            read_coord_names = fieldnames(dcm.coords);
            for kk = 1:numel(cntrl_coord_names)
                chk_nm = cntrl_coord_names{kk};
                diff = diff + ...
                    double(~any(strcmpi(chk_nm,read_coord_names)));
            end 

            testCase.assertLessThan(diff, 5e-4,"Disagreement b/t " + ...
                "quaternion-calc'd & HPR-calc'd orientmat");
        end

    end

    methods (Static)

        function diff = check_hpr(h, p, r, omatin)
            omat = euler2orient(h, p, r);
            diff = sum(omat - omatin,"all","omitnan");
            [heading, pitch, roll] = orient2euler(omat);
            diff = diff + sum(heading - h,"all","omitnan");
            diff = diff + sum(pitch - p,"all","omitnan");
            diff = diff + sum(roll - r,"all","omitnan");
        end

    end

end

