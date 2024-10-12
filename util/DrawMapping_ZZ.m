function DrawMapping_ZZ(monkey_hemi)
% DrawMapping.m
% First version by HH @ Gu Lab 2013

monkey_hemis = {'Fara_L','Fara_R','Messi_R','Messi_L','Duncan_L','Duncan_R'};

if nargin == 0
        monkey_hemi = 1;
%     monkey_hemi = 2;
end
clc; clear global

%% Parameters
global maxX maxY maxZ movDis GMTypes data;
global hemisphere monkey; global gridRange;
global MRI_path MRI_offset AP0; global MRI_offset_new;
global linWid start_end_markers overlapping Duncan_left_AP0 Messi_left_AP0 Fara_left_AP0;
global overlapTransparent;

linWid = 1.3;
overlapping = [-15 15; -15 15]; start_end_markers = 0; % First row for area annotation; second for unit annotation
% overlapping = [0 0; 0 0]; start_end_markers = false; % First row for area annotation; second for unit annotation
% overlapping = [-2 2; -2 2]; start_end_markers =0; % First row for area annotation; second for unit annotation
overlapTransparent = 0.3;

maxX = 65; % Grid size
maxY = 35;
maxZ = 400;   % Drawing range (in 100 um)
movDis = 150 * 2; % Maximum travel distance of the microdriver (in 100 um)

% V = ones(65,35,maxZ + 50,3);  % Matrix to be rendered. The z-resolution = 100 um.
V = ones(maxX,maxY,maxZ + 50,3);  % Matrix to be rendered. The z-resolution = 100 um.
Duncan_left_AP0 = 56; % For MRI alignment. This controls MRI-AP alignment so will affect all monkeys, whereas "AP0" of each monkey below only affects AP-Grid alignment.  HH20150918
Duncan_left_AC0 = 31;
% row 30th is AC-1 (positive indicates anterior, while negtive is posterior)

Messi_left_AP0 = 40;
Messi_left_AC0 = 17;

Fara_left_AP0 = 47;
Fara_left_AC0 = 19;
%% Define our area types here
GMTypes = { % 'Type', ColorCode (RGB);
    'GM',[0.3 0.3 0.3];    % Gray matter without visual modulation, lightgray
    %     'VM',[0 0 0];          % Visual modulation but not sure to which area it belongs
    'CD',[1 0 0];    % caudate
    'GP',[1 1 1];        % globus pallidus and so on
    'Pu',[1 1 1];    % putamen
    %     'VIP',[1 0.65 0];
    'FC',[0.3 0.3 0.3];         % frontal cortex including F7, F6, 6DR(dorsal premotor area 6) and so on.
    'LV',[0 0 0];       % lateral ventricle
%     'CC',[0.1 0.1 0.1];  % cingulate cortex, darkgray
    'CC',[0 1 0];  % cingulate cortex, darkgray
    'NAc',[0 0 1];
    };
for i = 1:length(GMTypes)
    eval([GMTypes{i,1} '= -' num2str(i) ';']);  % This is a trick.
end

%% Input our mapping data here

switch monkey_hemis{monkey_hemi}
    
    case 'Fara_L'
        
%         toPlotTypes3D = [CD GP Pu FC LV CC];
%         toPlotTypes_zview = [ CD GP Pu LV ];
        toPlotTypes3D = [CD GM CC];
        toPlotTypes_zview = [ CD CC ];
        
        gridRange = [1 25; 1 20];
        monkey = 13;  hemisphere = 1;
        AP0 = Fara_left_AP0;
        
        data = {
            %{[Session(s)], [LocX(Posterior) LoxY(Lateral)], [GuideTube(cm) Offset(cm)], [AreaType, Begin(100um), End(100um); ...] , vode retrieval}
            % When you are not sure about one area, use "AreaType-100" instead
            {272, [15, 8], [2.3, 0.0], [GM 8 35; CD 141  165]}  % only to 165
            {273, [15, 6], [2.3, 0.0], [GM 0 26; CC 68 97; CD 142 170]} % only to 165
            {274, [15, 10], [2.4, 0.0], [GM 13 35; CD 132 145]}   % only to 143
            {275, [15, 9], [2.4, 0.0], [GM 8 35; CD 134 150]}    % only to 150
            {276, [15, 7], [2.3, 0.0], [GM 15 35; CD 137 160]} % only to 160
            {278, [17, 8], [2.3, 0.0], [GM 10 30; CD 145 160]} % only to 160
            {279, [19, 8], [2.3, 0.0], [GM 8 30; CD 142 180]} % only to 175
            {280, [16, 8], [2.3, 0.0], [GM 11 39; CD 147 180]} % only to 176
            {281, [16, 7], [2.4, 0.0], [GM 14 32; CD 138 170]} % only to 163
            {282, [16, 9], [2.4, 0.0], [CD 145 170]} % only to 165
            {283, [16, 6], [2.3, 0.0], [GM 12 36; CC 78 97; CD 154 185]} % only to 181
            {284, [17, 9], [2.3, 0.0], [GM 22 43; CD 159 195]} % only to 190
            {285, [18, 7], [2.3, 0.0], [GM 20 39; CD 158 175]} % only to 173
            {286, [18, 8], [2.3, 0.0], [GM 80 92; CD 150 165]} % only to 160
            {287, [18, 6], [2.3, 0.0], [CC 67 95; CD 146 170]} % only to 169
            {288, [18, 9], [2.3, 0.0], [GM 12 38; CD 145 165]} % only to 162
            {289, [19, 9], [2.3, 0.0], [GM 13 31; GM 67 100; CD 148 160]} % only to 158
            {290, [19, 7], [2.4, 0.0], [GM 0 20; GM  57 83; CD 138 165]} % only to 160
            {291, [20, 7], [2.3, 0.0], [GM 15 31; CD 148 170]} % only to 168
            {292, [20, 8], [2.3, 0.0], [GM 13 42; CD 145 165]} % only to 160
            {293, [20, 6], [2.3, 0.0], [GM 7 20; CC 64 96; CD 150 182; GM 211 217]} 
            {294, [21, 8], [2.3, 0.0], [GM 13 35; CC 78 93; CD 156 180]} % only to 175
            {295, [21, 7], [2.3, 0.0], [GM 11 27; CC 66 102; CD 159 186]} 
            {296, [21, 9], [2.3, 0.0], [GM 16 39; CD 146 165]}   % only to 160
            {297, [21, 10], [2.3, 0.0], [GM 19 47; CD 151 170]} 
            {298, [19, 10], [2.3, 0.0], [GM 16 39; CD 146 175]} 
            {299, [22, 7], [2.3, 0.0], [GM 23 31; GM 73 107; CD 157 190]} 
            {300, [22, 8], [2.3, 0.0], [GM 15 33; CD 151 190]}  
            {301, [17, 7], [2.3, 0.0], [CC 81 97; CD 156 170]}  % only to 168
            {302, [17, 10], [2.4, 0.0], [GM 81 97; CD 155 171]}  
            {303, [17, 6], [2.4, 0.0], [CC 81 96; CD 148 180]}  % only to 176
            {304, [14, 8], [2.4, 0.0], [GM 33 52; CD 156 175]}  % only to 172
            {305, [14, 7], [2.4, 0.0], [GM 17 45; CD 154 175]}  % only to 173
            {306, [14, 6], [2.4, 0.0], [GM 27 49; CD 152 160]}  % only to 156
            {307, [14, 9], [2.4, 0.0], [GM 16 32; CD 146 155]}  % only to 152
            {308, [13, 8], [2.4, 0.0], [GM 24 40; CD 150 165]}  % only to 163
            {309, [13, 9], [2.4, 0.0], [GM 22 41; CD 152 165]}  % only to 161
            {310, [13, 7], [2.4, 0.0], [GM 18 38; CD 147 160]}  % only to 159
            {311, [13, 6], [2.4, 0.0], [GM 11 36; CD 147 165]}  % only to 161
            {312, [12, 8], [2.4, 0.0], [GM 12 36; CD 149 165]}  % only to 160
            {313, [12, 6], [2.4, 0.0], [GM 17 42; CD 146 150]}  % only to 148
            {314, [10, 7], [2.4, 0.0], [GM 12 37; CD 149 170]}  % only to 168
            {315, [11, 7], [2.4, 0.0], [GM 13 41; CD 148 160]}  % only to 154
            {316, [12, 7], [2.4, 0.0], [GM 15 40; CD 153 170]}  % only to 169
            {317, [12, 5], [2.4, 0.0], [GM 12 47; CD 156 175]}  % only to 172
            {318, [11, 8], [2.4, 0.0], [GM 25 48; CD 157 180]}  % only to 177
            {319, [11, 6], [2.4, 0.0], [GM 14 42; CD 146 155]}  % only to 153
            {320, [11, 9], [2.4, 0.0], [GM 17 41; GM 79 92; CD 151 160]}  % only to 156
            {321, [10, 6], [2.4, 0.0], [GM 17 37; CD 149 180]}  % only to 175
            {322, [10, 8], [2.4, 0.0], [GM 22 37; CD 154 175]}  % only to 172
            {323, [10, 5], [2.4, 0.0], [GM 20 41; CC 75 102; CD 151 180]} % only to 175
            {324, [9, 8], [2.4, 0.0], [GM 22 39; GM 89 107; CD 155 180]}  % only to 178
            }';
        
        MRI_path = 'Z:\Data\MOOG\Fara\Mapping\MRI\FaraOutput\forDrawMapping\';
        MRI_offset = {[-70.9547, 67.0547], [-233.02, 565.02],[0, -2.5]}; % Manual adjust HH20180604
        
        MRI_offset_new = {[-3.45, 191, 245.083],[0.087, -0.004]};   % Grid of Fara has 5 degree slope
        MRI_offset_new = {[-3.45, 196, 245.083],[0.087, -0.004]};
        
    case 'Fara_R'
        
%         toPlotTypes3D = [CD GP Pu FC LV CC];
%         toPlotTypes_zview = [ CD GP Pu LV ];
        toPlotTypes3D = [CD GM CC];
        toPlotTypes_zview = [ CD CC ];

        
        gridRange = [1 25; 1 20];
        monkey = 13;  hemisphere = 2;
        AP0 = Fara_left_AP0;
        
        data = {
            %{[Session(s)], [LocX(Posterior) LoxY(Lateral)], [GuideTube(cm) Offset(cm)], [AreaType, Begin(100um), End(100um); ...] , vode retrieval}
            % When you are not sure about one area, use "AreaType-100" instead
            {258, [10, 6], [1.9, 0.0], [GM 40 57; GM 123 143; GM 180 193]}  % only to 193
            {261, [15, 7], [2.3, 0.0], [GM 6 18; GM 71 80; CD 133 190]}
            {265, [15, 9], [2.3, 0.0], [GM 0 23; CD 130 166]}
            {267, [15, 10], [2.3, 0.0], [GM 8 26; GM 69 87; CD 133 167]}
            {268, [17, 7], [2.3, 0.0], [GM 8 27; CD 144 172]} % only to 172
            {269, [10, 8], [2.3, 0.0], [GM 6 24; CD 135 172]} % only to 171
            {270, [13, 8], [2.3, 0.0], [GM 5 12; CD 136 165]} % only to 164
            {271, [19, 7], [2.3, 0.0], [GM 57 87; CD 145 165]} % only to 163
            }';
        
        MRI_path = 'Z:\Data\MOOG\Fara\Mapping\MRI\FaraOutput\forDrawMapping\';
        MRI_offset = {[-70.9547, 67.0547], [-233.02, 565.02],[0, -2.5]}; % Manual adjust HH20180604
        
        MRI_offset_new = {[-0.95, 191, 245.083],[0.087, -0.004]};     % Grid of Fara has 5 degree slope
        
        
    case 'Messi_R'
        
%         toPlotTypes3D = [CD GP Pu FC LV CC];
%         toPlotTypes_zview = [ CD GP Pu LV ];
        toPlotTypes3D = [CD GM CC];
        toPlotTypes_zview = [ CD CC ];

        
        gridRange = [1 30; 1 20];
        monkey = 10;  hemisphere = 2;
        AP0 = Messi_left_AP0;
        
        data = {
            %{[Session(s)], [LocX(Posterior) LoxY(Lateral)], [GuideTube(cm) Offset(cm)], [AreaType, Begin(100um), End(100um); ...] , vode retrieval}
            % When you are not sure about one area, use "AreaType-100" instead
            {234, [17, 5], [2.3, 0.0], [GM 45 58; GM 70 84; CD 128 150]}  % only to 150
            {235, [17, 8], [2.3, 0.0], [GM 14 30; GM 66 93; CD 134 183; LV 189 199]}
            {236, [17, 6], [2.8, 0.0], [GM 0 26; GM 39 47; CD 77 94]} % only to 94
            {237, [18, 6], [2.7, 0.0], [GM 34 62; CD 88 106]} % only to 90    % 106 is gotten from session 238
            {239, [19, 6], [2.7, 0.0], [GM 11 64; CD 91 112]} % only to 112
            {240, [18, 5], [2.7, 0.0], [GM 10 61]} % only to 61
            {241, [18, 5], [2.7, 0.0], [CD 82 120]} % only to 120
            {242, [17, 7], [2.7, 0.0], [CD 80 90]} % only to 85
            {243, [17, 7], [2.6, -0.1], [GM 36 48; GM 59 67; CD 103 128]} % only to 128
            {244, [18, 7], [2.6, -0.1], [GM 48 73; CD 115 148]}
            {245, [18, 8], [2.6, -0.1], [GM 51 66; CD 105 132]} % only to 131
            {246, [18, 9], [2.6, -0.1], [GM 7 27; CD 104 140]}
            {247, [19, 8], [2.6, -0.1], [GM 49 71; CD 107 155]} % only to 126    % 155 is gotten from session 248
            {249, [19, 9], [2.6, -0.1], [GM 4 24; CD 106 140]}   % only to 114     % 140 is gotten from session 250
            {251, [15, 8], [2.6, -0.1], [GM 5 21; CD 114 130]}   % only to 130
            {252, [15, 9], [2.6, -0.1], [GM 15 34; CD 114 135]} % only to 135
            {254, [15, 7], [2.8, -0.05], [GM 19 59; CD 83 100]} % only to 98
            {255, [16, 8], [2.8, -0.1], [GM 25 49; CD 85 100]}  % only to 98
            {256, [16, 9], [2.8, -0.1], [GM 22 41; CD 83 110]}
            }';
        
        MRI_path = 'Z:\Data\MOOG\Messi\Mapping\MRI\MessiOutput\forDrawMapping\';
        MRI_offset = {[-70.9547, 67.0547], [-233.02, 565.02],[0, -2.5]}; % Manual adjust HH20180604
        
        MRI_offset_new = {[-1.45, 321, 245.083],[0, 0]};
        MRI_offset_new = {[-1.45, 351, 245.083],[0, 0]};   % align to MRI on 20221208
        % MRI_offset_new = {[1.55, 191, 254.984],[0, -2.5]};   % after manually align with DrawMapping_ZZ
        
        
    case 'Messi_L'
        
%         toPlotTypes3D = [CD GP Pu FC LV CC];
%         toPlotTypes_zview = [ CD GP Pu LV ];
        toPlotTypes3D = [CD GM CC];
        toPlotTypes_zview = [ CD CC ];

        
        gridRange = [1 30; 1 20];
        monkey = 10;  hemisphere = 1;
        AP0 = Messi_left_AP0;
        
        data = {
            %{[Session(s)], [LocX(Posterior) LoxY(Lateral)], [GuideTube(cm) Offset(cm)], [AreaType, Begin(100um), End(100um); ...] , vode retrieval}
            % When you are not sure about one area, use "AreaType-100" instead
            {216, [17, 8], [2.2, 0.0], [GM 22 90; CD 130 161]}  % only to 161
            {217, [17, 8], [2.4,-0.05], [GM 20 85; CD 106 120]}
            {218, [17, 4], [2.3,0], [GM 16 74; CD 127 135; GM 152 167; GM 185 199]} % only to 199
            {219, [19, 7], [2.4, 0], [GM 23 72; CD 128 148]}  % only to 148
            {220, [16, 5], [2.3, 0], [GM 13 34; GM 60 68; GM 141 150]}  % only to 150
            {221, [17, 6], [2.3, 0], [GM 6 63; CD 123 152; GM 173 186; GM 207 220]}
            {222, [17, 5], [2.3, 0], [GM 13 19; GM 31 65; CD 133 142]} % only to 142
            {223, [18, 7], [2.3, 0], [GM 21 41; GM 58 82; CD 135 153]} % only to 153
            {224, [18, 8], [2.4, 0], [GM 38 74]}   % non-ideal signals in deeper position
            {225, [18, 6], [2.4, 0], [GM 24 37; GM 64 79; CD 139 147]} % only to 147, still non-ideal signals in deeper position
            {226, [19, 6], [2.3, 0], [GM 58 92; CD 123 162; GM 169 181]} % to 198
            {228, [19, 8], [2.8, 0], [GM 0 17; CD 72 87; GM 104 120; GM 129 141]}
            {229, [20, 7], [2.7, 0], [GM 5 46; CD 55 63]} % only to 63
            {230, [16, 8], [2.9, 0], [GM 6 50; CD 58 66]} % only to 66
            {231, [15, 8], [2.9, 0], [GM 4 48; CD 54 10; GM 115 131]} % only to 131
            {232, [15, 7], [2.9, 0], [GM 2 46; CD 50 61; CD 109 110]} % only to 110
            }';
        
        MRI_path = 'Z:\Data\MOOG\Messi\Mapping\MRI\MessiOutput\forDrawMapping\';
        MRI_offset = {[-70.9547, 67.0547], [-233.02, 565.02],[0, -2.5]}; % Manual adjust HH20180604
        
        MRI_offset_new = {[-1.45, 321, 245.083],[0, 0]};
        % MRI_offset_new = {[1.55, 191, 254.984],[0, -2.5]};   % after manually align with DrawMapping_ZZ
        
        
    case 'Duncan_L'
        
        %  Duncan_Left
        %  %{
        
        % Header
%         toPlotTypes3D = [ CD GP Pu FC LV CC GM NAc];    % Which area types do we want to plot in 3D plot?
%         %         toPlotTypes3D = [ CD GP Pu VIP FC LV CC];    % Which area types do we want to plot in 3D plot?
%         toPlotTypes_zview = [ CD GP Pu LV GM NAc];    % Which area types do we want to plot ?
%         %         toPlotTypes_zview = [ CD GP VIP Pu LV ];    % Which area types do we want to plot ?
        
        toPlotTypes3D = [CD GM CC];
        toPlotTypes_zview = [ CD CC ];

        
        gridRange = [15 45; 1 20];  % Grid plotting ragne = [xLow xHigh; yLow yHigh]
        
        monkey = 15;
        hemisphere = 1;
        AP0 = Duncan_left_AP0 ;
        
        data = {
            %{[Session(s)], [LocX(Posterior) LoxY(Lateral)], [GuideTube(cm) Offset(cm)], [AreaType, Begin(100um), End(100um); ...] , vode retrieval}
            % When you are not sure about one area, use "AreaType-100" instead
            %             {1,[46,17],[1.8,0.0],[GM 0 100]}
            %             {1,[55,11],[1.8,0.0],[GM 0 210]}
            %
            %             {1,[46,15],[1.7,0.0],[FC 33 81]}
            %             {1,[46,13],[1.7,0.0],[FC 24 93; GM 182 210]}
            %             {2,[53,15],[1.9,0.0],[FC 21 65; VIP 84 113;GM 229 284]}
            %             {3,[53,13],[1.8,0.0],[FC 27 43; CC 67 92; VIP 110 130; GM 246 260]}
            %             {4,[55,13],[1.8,0.0],[FC 24 50; VIP 105 132; GM 242 261]}
            %             {5,[55,11],[1.8,-0.05],[FC 25 45; CC 65 102; VIP 120 150; GM 212 268]}
            %             {6,[55,9],[1.8,0.0],[FC 21 36; VIP 90 123]}
            %             {10,[30,9],[1.4,0.0],[FC 35 73; CD 171 220]}
            %             {11,[30,7],[1.4,-0.2],[FC 60 81; CC 115 141; CD 202 231]}    % 从guidetube内距离尖端2mm处出发
            %             {12,[30,11],[1.8,0.0],[FC 21 42; CD 185 216; GP 224 265; GM 280 290]}
            %手动测量得到现在的微操前1000um没动，前2000um少走~1200um，前3000um少走~1500um，4000um到14000um走了~10500um。所以最终下15000um总共少走~1000um
            %             {10,[30,9],[1.4,0.0],[FC 20 48; CD 148 195]}
            %             {11,[30,7],[1.4,-0.2],[FC 45 66; CC 103 131; CD 177 206]}    % 从guidetube内距离尖端2mm处出发
            %             {12,[30,11],[1.8,0.0],[FC 6 27; CD 160 191; GP 200 240; GM 260 270]}
            {13, [30,5], [1.75,-0.1], [FC 23 36; CC 58 80.5; CC 86 104; CD 138 195; GP 210 230]}  %约210进入GP 但是不知道多少出230，这里为了不报错随便写了一个值
            {15, [30,7], [1.5, 0.0], [FC 14 36; CC 76 101; CD 152 190; GM 221 241]} % till 241
            {16, [30,11],[1.75 -0.1], [FC 21 48; CD 150 163; Pu 195 207; GP 219 239]}
            {17, [30,6], [1.5 0.0], [FC 15 40; CC 73 97; CC 92 108; CD 149 197; GP 230 250]}
            {19, [30,8], [1.6,0.0], [FC 9 34; CD 136 167]} % 只走到了167
            {20, [29,7], [1.6,0.0], [CC 76 105; CD 144 172]} % 只到了172
            {21, [32,7], [1.75,0.0],[FC 9 23; CC 61 91; CD 127 149]} % only to 149
            {22, [35,8], [1.75,0.0],[FC 4 30; CC 66 83; CD 135 158]} % 158 out of CD
            {23, [38,9], [1.75,0.0],[FC 7 30; CD 128 141]}
            {24, [40,9], [1.75,0.0],[FC 7 29; CD 136 148]}
            {25, [30,9], [1.75,0.0],[FC 10 28; CD 119 150]}  % only to 144
            {26, [28,6], [1.75,-0.1],[FC 13 31; CC 57 90; CD 128 170; GP 209 229]}
            {27, [26,6], [1.6,0.0], [FC 10 31; CC 63.5 89; CD 142 168]}    % only to 168
            {28, [22,6], [1.6,0.0], [FC 24 48; CC 73.7 111; CD 156 220]}   % CD= CD+NAc only to 220
            {29, [20,6], [1.75,0.0],[FC 3 25; CC 52 84; CD 132 193]} % only to 193
            {30, [31,6], [2.0,0.0], [CC 10 28; CC 33.5 48; CD 93 109]} % only to 109
            {31, [31,7], [1.75,0.0],[FC 3 31; CC 47 64.5; CC 70 86; CD 131 168]} % only to 168
            {32, [31,8], [1.75,0.0],[CC 78.5 88; CD 136 178]}
            {33, [31,9], [1.6, 0.0],[FC 16 38; CD 139 155]} % only to 155
            {34, [30,10],[1.6, 0.0],[FC 21 51; CD 140 166]}
            {35, [20,7], [1.6, 0.0],[FC 15 43; CD 151 164]} % only to 164
            {36, [20,8], [1.6, 0.0],[FC 21 57; CD 150 160]} % only to 160
            {37, [20,9], [1.6, 0.0],[FC 15 44; CD 141 176]} % only to 176
            {38, [20,10],[1.6, 0.0],[FC 13 49; FC 78 93; CD 149 167]} % only to 167
            {39, [19,7], [1.6, 0.0],[FC 14 50; CC 78 101;LV 149 154; CD 155 181]} % only to 181
            {40, [19,8], [1.6, 0.0],[FC 11 45; CD 147 167]} % only to 167
            {41, [19,9], [1.6, 0.0],[FC 18 45; FC 96 97; CD 145 179]} % only to 179
            {42, [19,10],[1.6, 0.0],[FC 21 48; CD 146 185; GM 208 228]} % only to 208
            {43, [18,7], [1.6, 0.0],[FC 11 42; CC 79 95; CD 157 205; GM 217 237]} % only to 225
            {44, [18,8], [1.6, 0.0],[FC 10 34; GM 63 95; LV 142 155; CD 159 200]} % only to 200
            {45, [18,9], [1.6, 0.0],[FC 20 63; CD 167 196]} % only to 196
            {46, [18,10],[1.6, 0.0],[FC 9 40; GM 47 71; CD 145 169]} % only to 169
            {47, [17,7], [1.75,0.0],[FC 3 25; CC 58 82; LV 135 140; CD 140 152]} % only to 152
            {48, [17,8], [1.75,0.0],[FC 7 29; CC 50 93; LV 143 160; CD 162 171; GM 211 223]}
            {49, [17,9], [1.75,0.0],[FC 3 26; CD 133 154]} % only to 154
            {50, [17,10],[1.75,0.0],[FC 7 33; CD 141 153]} % only to 153
            {51, [17,11],[1.8, 0.0],[FC 1 25; CD 129 135]} % only to 135
            {52, [17,6], [1.75,0.0],[FC 10 33; CC 53 73; CC 80 93; LV 142 158; CD 158 174; GM 197 199]} % only to 199
            {53, [18,6], [1.75,0.0],[FC 4 28; CC 58 86; LV 129 134; CD 135 150]} % only to 150
            {54, [17,12],[1.8, 0.0],[FC 5 77]}
            %             {54, [17,11],[1.7,0.0], [FC 22 53; GM 206 220]}
            %             {55, [15,6], [1.8,0.0], [FC 11 35; CC 55 72; CC 77 93; LV 141 150; CD 155 221]}
            %             {56, [15,7], [1.75,0.0], [FC 21 45; CC 66 99; CD 160 178; GM 192 212]}
            {57, [20,5], [1.75,0.0], [FC 14 35; CC 57 76; CC 83 94; CD 137 211]} % VM in 23 but not in near holes; CD = CD + NAc; only to 211
            {58, [20,4], [1.75,0.0], [FC 15 42; CC 60 81; CC 85 100; LV 137 148; CD 148 180; LV 183 190; NAc 190 231]} % no VM; only to 231; NAc ?
            {58, [20,11],[1.75,0.0], [FC 22 98; CD 145 160; Pu 176 196]} % only to 196
            {59, [20,3], [1.75,0.0], [FC 23 50; CC 63 87; CC 96 110; LV 149 155; CD 155 177; LV 177 194; NAc 194 212]} % only to 212
            {60, [19,5], [1.8, 0.0], [FC 23 43; CC 60 76; CC 84 98; LV 133 145; CD 146 187]} % CD = CD + NAc; only to 187
            {61, [19,11],[1.8, 0.0], [GM 42 65; GM 83 98; CD 147 174; Pu 187 206]} % only to 206
            {62, [19,6], [1.6, 0.0], [CC 68 83; CC 107 123; CD 181 201]} % only to 181
            {62, [20,7], [1.6, 0.0], [CC 56 75; CD 155 199]} % only to 199
            {63, [21,7], [1.6, 0.0], [CC 59 79; CC 101 117; CD 160 187]} % only to 187
            {64, [21,9], [1.6, 0.0], [CC 71 87; CD 160 176]} % only to 176
            {65, [21,8], [1.6, 0.0], [CC 75 90; CD 167 193]} % only to 193
            {66, [21,6], [1.6, 0.0], [CC 72 99; CC 108 129; CD 164 210]}
            {67, [21,10],[1.6, 0.0], [GM 81 101;CD 167 195]} % only to 195
            {68, [25,7], [1.6, 0.0], [CC 64 101; CD 163 191]} % only to 191
            {69, [25,9], [1.6, 0.0], [GM 74 98; CD 161 183]} % only to 183
            {70, [27,7], [1.6, 0.0], [CC 76 108; CD 186 211]} % only to 211   I'm not sure whether the mapping is right
            {71, [28,7], [1.6, -0.2], [CC 84 112; CD 188 204]} % only to 204
            {72, [25,8], [1.6, 0.0], [CC 84 110; CD 171 197]} % only to 197
            {73, [24,8], [1.6, 0.0], [CC 88 108; CD 167 185]} % only to 185
            {74, [26,7], [1.6, 0.0], [CC 83 110; CD 174 200]} % only to 200
            {75, [24,7], [1.6, -0.1],[CC 101 121; CD 180 195]} % only to 195
            {76, [27,8], [1.6, -0.1],[CC 82 108; CD 185 207]} % only to 207
            {77, [35,7], [1.7, 0.0], [CD 140 168]} % only to 168
            {78, [37,8], [1.7, -0.1],[CC 83 108; CD 150 180]} % only to 180
            {79, [23,7], [1.6, 0.0], [CC 56 78; CD 156 186]} % only to 186
            {80, [33,7], [1.8, 0.0], [FC 13 39; CC 70 104; CD 138 162]}
            {81, [36,8], [1.8, 0.0], [FC 18 48; CC 82 106; CD 147 177; GP 183 203]} % only to 183
            {82, [36,7], [1.8, 0.0], [FC 17 46; CC 77 108; CD 147 174; GP 184 204]} % only to 187
            {83, [36,6], [1.8, 0.0], [FC 13 40; CC 69 87; CC 94 111; CD 148 173; GP 179 228]} % till 228
            {84, [36,9], [1.9, 0.0], [FC 9 36; CD 139 151]} %only to 151
            {85, [35,7], [1.8, 0.0], [FC 16 43; CC 75 90; CC 96 110; CD 148 180; GP 186 206]}
            {86, [32,8], [2.05,0.0], [FC 8 26; CD 105 117]} % only to 117
            {87, [32,9], [2.15,0.0], [CD 103 115]}  %only to 115
            {88, [32,10],[2.05,0.0], [CD 109 133]} % only to 133
            {89, [33,8], [2.15,0.0], [CD 95 127]} % only to 127
            {90, [33,10],[2.05,0.0], [CD 109 124]} % only to 124
            {91, [35,10],[2.05,0.0], [CD 113 127]} % only to 127
            {92, [36,10],[2.05,0.0], [CD 115 140; Pu 149 169]} % 149 is around the next GM
            {93, [34,8], [2.05,0.0], [CD 113 135; Pu 143 163]} % only to 148
            {94, [34,10],[2.05,0.0], [CD 110 128]} % only to 128
            {95, [22,9], [2.15,0.0], [CD 103 132]} % only to 132
            {96, [22,10],[2.05,0.0], [CD 119 145]} % only to 145
            {97, [22, 8], [2.15, 0.0], [CD 102 142]} % only to 142
            {98, [23,11],[2.05, 0.0], [CD 117 142]} % only to 142
            {99, [24,10],[2.05, 0.0], [CD 104 145; GM 171 180]} % only to 177
            {100, [24,11],[2.05,0.0],[GM 59 79; CD 116 121]} % only to 121
            {101, [25,11],[2.05,0.0],[CD 112 122]} % only to 112
            {102, [26,11],[2.05,0.0],[CD 112 142; GM 159 169]} % only to 166
            {103,[26, 10],[2.05,0.0],[GM 12 33; GM 65 88; CD 121 135]} % only to 135
            {104,[27, 11],[2.05,0.0],[GM 5 33; CD 113 132]} % only to 132
            {105,[28, 11],[2.05,0.0],[GM 11 40; CD 114 138; GM 151 161]} % only to 152
            {117,[55, 16],[2.05,0.0],[GM 16 40; GM 70 90]}  % only to 90, modest VM in 7988
            {118,[53, 18],[2.05,0.0],[GM 23 49; GM 57 87]} % 10382 colse to MST
            {120,[55, 18],[2.15,0.0],[GM 10 42]}
            {123, [53, 18],[2.25,0.0],[GM 38 56]}  % only to 5080
            {124, [53, 19],[2.25,0.0], [GM 17 43]}  % only to 4227
            {125, [53, 19],[2.25,0.0], [GM 17 57; GM 91 95]}  %only to 9514
            {126, [53, 16],[2.15,0.0], [GM 14 52; GM 71 86]}
            {127, [28, 10],[2.3,0.0], [CD 62 95]}       %%%%%%  Microstimulation From here
            {129, [32, 9],[2.0 0.0], [GM 0 16; CD 93 115]}
            {130, [30,8], [1.8 0.0], [GM 11 20; CD 119 132]}
            {131, [17,9], [1.8,0.0], [GM 2 23; CD 129 150]}
            {132, [20,10], [2.0,0.0],[CD 101 125]} % only to 12243
            {133, [18,10],[2.0,0.0],[GM 0 22; CD 105 118]} %only to11719
            {134, [19,10],[2.0,0.0], [CD 100 120]} % only to 11548
            {135, [34,10],[2.0,0.0], [GM 2 35; CD 94 110]} % only to 10948
            {136, [34,10],[2.0,0.0], [CD 94 110]}  % only to 10780
            {137, [35,10],[1.8,0.0],[CD 128 142]} % only to 14187
            {138, [30,10],[2.0,0.0],[GM 14 21; CD 106 124]} % only to 12366
            {139, [25,9],[2.0,0.0],[CD 102 130]}
            {140, [24,10],[2.0,0.0],[GM 0 24; CD 99 120]} % only to 11942
            {141, [26,10],[2.0,0.0],[GM 6 15; CD 104 121]} % only to 12067
            {142, [29,10],[2.0,0.0],[CD 97 121]}
            {143, [29, 8],[2.0,0.0],[CC 25 50; CD 86 100]} % only to 9723
            {144, [27, 8],[2.0,0.0],[CC 33 56; CD 97 110]} % only to 10882
            {145, [28,7], [1.8, 0.0], [CC 58 81;CD 122 135]} % only to 13346
            {146, [25, 8],[2.0,0.0], [CD 104 132]}
            {147, [28,7],[2.0,0.0],[CC 27 54; CD 94 105]} % only to 10352
            {148, [27,9],[2.0,0.0],[CD 91 110]} % only to 10763
            {149, [26,10],[2.0,0.0],[CD 110 138]}
            {150, [24,8],[1.8,0.0],[GM 0 23; CD 110 135]} % only to 13420
            {151, [25,7],[1.8,0.0],[GM 2 26; CD 117 135]} % only to 13384
            {152, [30,9],[1.8,0.0],[GM 0 23; CD 111 125]} % only to 12140
            {153, [22,10],[2.0,0.0],[GM 0 18; GM 42 60; CD 110 120]} % only to11488
            {154, [20,9],[2.0,0.0],[CD 98 120]} % only to 11760
            {155, [18,8],[1.9,0.0],[GM 1 19; CD 119 146 ]}
            {156, [17,10],[2.0,0.0],[CD 122 143]}
            {157, [30,7],[1.8,0.0],[CC 55 70; CD 124 150]}
            {158, [32,10],[1.9,0.0],[GM 0 19; CD 109 138]}
            {159, [31,10],[1.8,0.0],[CD 122 145]}
            {160, [26,8],[1.8,0.0],[CD 132 150]}
            {161, [22,9],[2.0,0.0],[CD 106 120]} %only to 11878
            {162, [28,9],[1.9,0.0],[CD 102 125]} % only to 12049
            {163, [28,9],[2.0,0.0],[CD 93 120]} % only to 11171
            {165, [31,9],[1.9,0.0],[CD 111 125]} % only to 125
            {166, [31,9],[1.9,0.0],[CD 111 135]}
            {167, [27,11],[1.9,0,0],[CD 110 130]}
            {168, [24,7],[1.9,0.0],[GM 6 25; CD 114 136]}
            {169, [31,9],[1.9,0.0],[GM 3 25; CD 111 140]}
            {170, [32,8],[1.8,0.0],[GM 0 16; CC 57 73; CD 115 145]}
            {171, [29,9],[1.9,0.0],[GM 5 23; CD 108 125]} % only to 12329
            {172, [26,9],[1.9,0.0],[GM 0 19; CD 105 120]} % only to 11510
            {173, [23,10],[2.0,0.0],[GM 0 18; GM 48 62; CD 100 135]}
            {174, [20,8],[2.0,0.0],[CD 106 130]}
            {175, [24,9],[2.0,0.0],[CD 104 125]}
            {176, [25,10],[2.0,0.0],[CD 102 110]} % only to 10550
            {177, [25,10],[2.0,0.0],[CD 102 120]} % only to 11553
            {178, [27,10],[2.0,0.0],[CD 101 110]}
            {179, [27,10],[1.9,0.0],[GM 7 28; CD 113 125 ]} % only to 12249
            {180, [27,10],[1.9,0.1],[CD 101 120]}
            {181, [31,8],[1.9,0.0],[CD 120 125]} % only to 12473
            {182, [31,8],[1.9,0.0],[CC 44 71; CD 120 150]}
            {183, [33,10],[1.9,0.0],[GM 2 25; CD 110 125]} % only to 12085
            {184, [33,9],[1.9,0.0],[CD 113 120]} % only to 11803
            {185, [33,9],[1.9,0.0],[GM 7 17; CD 113 135]}
            {186, [34,9],[1.9,0.0],[GM 0 22; CD 114 125]} % only to 12397
            {187, [31,7],[1.8,0.0],[CC 58 87; CD 128 150]}
            {188, [33,7],[1.9,0.0],[GM 6 30; CC 43 61; CC 65 82; CD 121 150]}
            {189, [33,8],[1.9,0.0],[GM 3 15; CC 55 73; CD 115 140]}
            {190, [28,8],[1.9,0.0],[GM 7 21; CD 116 140]}
            {191, [21,9],[1.9,0.0],[GM 7 26; CD 115 130]}; %only to 12645
            {192, [21,10],[1.9,0.0],[GM 8 35; CD 122 140]}
            
            }';
        
        MRI_path = 'Z:\Data\MOOG\Duncan\Mapping\MRI\DuncanOutput\forDrawMapping\';
        %         MRI_offset = {[-68 79]*1.1-9.5 ,[-240 500]*1.1+23 , [0 -2.5]};   % [x1 x2],[y1 y2], [dx/dxSelect slope, dy/dxSelect slope]
        
        MRI_offset = {[-70.9547, 67.0547], [-233.02, 565.02],[0, -2.5]}; % Manual adjust HH20180604
        
        % Update version with ratio fixed {[x_center, y_center, x_range], [slopes]}. HH20180606
        %         MRI_offset_new = {[-1.9500,  166,  138.0094 ], [0, -2.5]};
        %         MRI_offset_new = {[-0.45, 236, 245.083],[0, -2.5]};
        %         MRI_offset_new = {[-0.95, 106, 254.984],[0, 0]};  % align with electrode
        %         MRI_offset_new = {[0.55, 116, 254.984],[0, 0]};   % align with the row 30
        %         MRI_offset_new = {[0.55, 126, 265.285],[0, 0]}; % align with the row 20
        MRI_offset_new = {[0.55, 116, 265.285],[0, 0]};  % align with row 31 which may be AC0
        
        
    case 'Duncan_R'
        
%         toPlotTypes3D = [ CD GP Pu FC LV CC];
%         toPlotTypes_zview = [ CD GP Pu LV ];
        toPlotTypes3D = [CD GM CC];
        toPlotTypes_zview = [ CD CC ];

        
        gridRange = [15 45; 1 20];
        
        monkey = 15; hemisphere = 2;
        AP0 = Duncan_left_AP0;
        
        
        data = {
            %{[Session(s)], [LocX(Posterior) LoxY(Lateral)], [GuideTube(cm) Offset(cm)], [AreaType, Begin(100um), End(100um); ...] , vode retrieval}
            % When you are not sure about one area, use "AreaType-100" instead
            {125, [30, 7], [1.7, -0.1], [FC 17 25; CC 55 77; CC 83 101; CD 133 187; GP 208 228]}  % the largest depth of GP is chosen artificially
            %             {14, [30,7], [1.7,-0.1], [FC 14 36; CC 76 101; CD 144 184; GP 210 230]}     % Right hemisphere
            }';
        
        MRI_path = 'Z:\Data\MOOG\Duncan\Mapping\MRI\DuncanOutput\forDrawMapping\';
        MRI_offset = {[-70.9547, 67.0547], [-233.02, 565.02],[0, -2.5]}; % Manual adjust HH20180604
        
        MRI_offset_new = {[-0.45, 236, 245.083],[0, 0]};
        % MRI_offset_new = {[1.55, 191, 254.984],[0, -2.5]};   % after manually align with DrawMapping_ZZ
        
end


%% ========== 3-D Visualization ========== %%

for channel = 1:length(data)
    gridLoc = data{channel}{2};
    GMData = data{channel}{4};
    if isempty(GMData); continue; end
    
    GuideTubeAndOffset = data{channel}{3};
    offSet = round((GuideTubeAndOffset(1) + GuideTubeAndOffset(2) - 1.4) * 100);  % Related to Guide Tube 1.4 cm!!
    %      offSet = (data{channel}{3}(2)+ data{channel}{3}(1))*100;
    
    
    for GMType = -size(GMTypes,1):-1  % For each area type
        GMThisType = find(GMData(:,1) == GMType);   % Read out ranges for this type
        if isempty(GMThisType) || isempty(intersect(toPlotTypes3D,GMType)); continue; end       % If absent in this channel or absent in areas we want to plot, next type
        
        GMPos = [];    % Clear cache
        for i = 1:length(GMThisType)    % For each appearance
            % Attach each appearance to position cache
            GMPos = [GMPos round(maxZ - offSet - GMData(GMThisType(i),3)):round(maxZ - offSet - GMData(GMThisType(i),2))];
            %              GMPos = [GMPos round(offSet + GMData(GMThisType(i),2)):round(offSet + GMData(GMThisType(i),3))];
            
        end
        
        % Add color to positions in the cache
        try
            if hemisphere == 2
                V(gridLoc(1), maxY - gridLoc(2), GMPos, :) = repmat(GMTypes{-GMType,2},length(GMPos),1);  %记右脑的时候再改
            else
                V(gridLoc(1), gridLoc(2), GMPos, :) = repmat(GMTypes{-GMType,2},length(GMPos),1);
            end
        catch
            fprintf('Warning: Out of Range (Session %g, Channel [%g,%g], GMType %s)\n', channel, data{channel}{2}(1), data{channel}{2}(2), GMTypes{-GMType,1});
            keyboard;
        end
    end
    
    % Add start and end markers
    %     if hemisphere == 2
    %         V(gridLoc(1), (maxY - gridLoc(2)),(maxZ - offSet - 2):(maxZ - offSet - 1),:) = repmat([0 0 0],2,1);
    %         if length(data{channel}) <5
    %             V(gridLoc(1), (maxY - gridLoc(2)),(maxZ - offSet - movDis - 1):(maxZ - offSet - movDis),:) = repmat([0 0 0],2,1);
    %         else
    %             V(gridLoc(1), (maxY - gridLoc(2)),(maxZ - offSet - data{channel}{5} - 1):(maxZ - offSet - data{channel}{5}),:) = repmat([0 0 0],2,1);
    %         end
    %     else
    %         V(gridLoc(1), (gridLoc(2)),(maxZ - offSet - 2):(maxZ - offSet - 1),:) = repmat([0 0 0],2,1);
    %         if length(data{channel})<5
    %             V(gridLoc(1), (gridLoc(2)),(maxZ - offSet - movDis - 1):(maxZ - offSet - movDis),:) = repmat([0 0 0],2,1);
    %         else
    %             V(gridLoc(1), (gridLoc(2)),(maxZ - offSet - data{channel}{5} - 1):(maxZ - offSet - data{channel}{5}),:) = repmat([0 0 0],2,1);
    %         end
    %     end
    
end

% Render the mapping results using "vol3d" function
% close all;
set(figure(801),'Position',[10 100 600 600]);  clf;
set(0, 'DefaultAxesXTickMode', 'auto', 'DefaultAxesYTickMode', 'auto', 'DefaultAxesZTickMode', 'auto');
vol3d('cdata',V);
view(-150,30);   % view(AZ,EL) sets the angle of the view from which an observer sees the current 3-D plot.

axis tight;
daspect([1,1,10]);  % sets the data aspect ratio.
alphamap('rampup'); %创建不透明度逐渐增加的线性 alphamap
alphamap(.5 .* alphamap);  % Transparency

grid on;  grid minor;
set(gca,'GridLineStyle','-');

set(gca,'xtick',-0.5:5:maxY-0.5);
if hemisphere == 2
    %     set(gca,'xticklabel','25|20|15|10|5|0');
    set(gca,'xticklabel',maxY:-5:0);
else
    %     set(gca,'xticklabel','0|5|10|15|20|25');
    set(gca,'xticklabel',0:5:maxY);
end
xlim([-1 maxY]);

set(gca,'ytick',-0.5:5:maxX-0.5);
% set(gca,'yticklabel','0|5|10|15|20|25|30');
set(gca,'yticklabel',0:5:maxX);
ylim([-1 maxX]);

set(gca,'ztick',0:50:maxZ+50);
% set(gca,'zticklabel','-250|-200|-150|-100|-50|0');
set(gca,'zticklabel',-maxZ-50:50:0);


ylabel('Posterior (X)')
xlabel('Lateral (Y)')

for i = 1:length(toPlotTypes3D)
    text(0, 0, maxZ-i*40, GMTypes{-toPlotTypes3D(i),1},'color',GMTypes{-toPlotTypes3D(i),2},'FontSize',20);
end

set(gcf,'color','w');
% set(findall(gcf,'fontsize',10),'fontsize',20);
SetFigure(20);

% k=1;
% for i = 1:150
%     view(-170+i*.5,20+i/10);
%     drawnow;
%     mov(k) = getframe(gcf);
%     k=k+1;
% end
%
% for i = 150:-1:1
%     view(-170+i*.5,20+i/10);
%     drawnow;
%     mov(k) = getframe(gcf);
%     k=k+1;
% end
%
% movie2avi(mov,'Test.avi');

%% ============ 2-D Visualization (Grid view from the top) =============== %

radius = 0.42;  % Radius of each hole (interval = 1)

% Plot grid outline
set(figure(802),'Position',[10 50 600 600]); clf
hold on; axis equal ij;  % ij: 将坐标轴设置为矩阵模式。此时水平坐标轴从左到右取值，垂直坐标从上到下

global h_grid_axis; h_grid_axis = gca;

x1 = gridRange(1,1); x2 = gridRange(1,2);
y1 = gridRange(2,1); y2 = gridRange(2,2);
xlim([y1-2,y2+2]);
ylim([x1-2,x2+2]);

% Quick select monkey_hemi. HH20160122
h_monkey_hemi = uicontrol('style','popupmenu','unit','norm','position',[0.01 0.94 0.131 0.035],...
    'string','Fara_L|Fara_R|Messi_R|Messi_L|Duncan_L|Duncan_R','callback',@change_monkey_hemi);
set(h_monkey_hemi,'value',monkey_hemi);

% Frame
interval = 5;
xLoc = intersect(x1:x2,0:interval:100);
yLoc = intersect(y1:y2,0:interval:100);
xLines = line(repmat([y1-1;y2+1],1,length(xLoc)),repmat(xLoc,2,1));  % xLines
yLines = line(repmat(yLoc+0.25,2,1), repmat([x1-1;x2+1],1,length(yLoc)));  % yLines
set(xLines,'LineWidth',5,'Color','g');
set(yLines,'LineWidth',5,'Color','g');
set(gca,'xtick', yLoc);
set(gca,'ytick', xLoc);

% Parallel drawing (to speed up)
xOffsets = repmat(x1:x2, y2-y1+1,1);
xOffsets = xOffsets(:)';
yOffsets = repmat([(y1:y2)+0.5*~mod(x1,2) (y1:y2)+0.5*mod(x1,2)], 1, fix((x2-x1 + 1)/2));
if mod(x2-x1+1,2) == 1
    yOffsets = [yOffsets y1:y2];
end
t = linspace(0,2*pi,100)';
xGrid = radius * repmat(sin(t),1,length(xOffsets)) + repmat(xOffsets,length(t),1);
yGrid = radius * repmat(cos(t),1,length(yOffsets)) + repmat(yOffsets,length(t),1);
set(fill(yGrid,xGrid,[1 1 1]),'LineWidth',1.5,'ButtonDownFcn',@SelectChannel);    % White color

% Plot mapping result
for channel = 1:length(data)
    xCenter = data{channel}{2}(1);
    yCenter =  data{channel}{2}(2) + 0.5 * ~mod(data{channel}{2}(1),2);
    
    % Paint
    haveTypes = intersect(toPlotTypes_zview,data{channel}{4}(:,1));    % Find which areas this channel has.
    if isempty(haveTypes)    % If there is no types of interest, paint using the color of GM
        t = linspace(0,2*pi,100)';
        xMap = radius * sin(t) + xCenter ;
        yMap = radius * cos(t) + yCenter;
        set(fill(yMap,xMap,GMTypes{-GM,2}),'ButtonDownFcn',@SelectChannel);
        %         c = 'k';
    else                     % Else, we paint colors for each type
        for iType = 1:length(haveTypes)
            t = linspace(2*pi/length(haveTypes)*(iType-1), 2*pi/length(haveTypes)*iType, 100)';
            xMap = 0.99 * radius * [0; sin(t)] + xCenter;
            yMap = 0.99 * radius * [0; cos(t)] + yCenter;
            set(fill(yMap,xMap,GMTypes{-haveTypes(iType),2}),'LineStyle','none','ButtonDownFcn',@SelectChannel);
        end
        
        %         if [0.21 0.72 0.07]*(GMTypes{-haveTypes(1),2})' > 0.5
        %             c = 'k';
        %         else
        %             c = 'y';
        %         end
    end
    
    %------ Denotations -----%
    % Session No.
    text(yCenter, xCenter-0.1, num2str(data{channel}{1}),'FontSize',7,'color','k','HorizontalAlignment','center','ButtonDownFcn',@SelectChannel);
    
end

% Legend
for i = 1:length(toPlotTypes_zview)
    text(y1+(i-1)*4, x1-3, GMTypes{-toPlotTypes_zview(i),1},'color',GMTypes{-toPlotTypes_zview(i),2},'FontSize',15);
end

set(gcf,'color','w')
set(gca,'fontsize',15)

if hemisphere == 1
    set(gca,'xdir','rev');
end

%% Load xls file
global num txt raw xls;
if isempty(num)
    XlsData = ReadXls('Z:\Labtools\ZZ_Tools\DataHub\DataHub.xlsm',2,3);
    num = XlsData.num;
    xls = XlsData.header;
    txt = XlsData.txt;
    
    disp('Xls Loaded');
end

%% Call back for 2-D Grid View
function SelectChannel(source,event)

global maxX maxY maxZ movDis GMTypes data Duncan_left_AP0 Messi_left_AP0 Fara_left_AP0
global MRI_path MRI_offset AP0; global MRI_offset_new;
global monkey hemisphere ;
global gridRange linWid start_end_markers overlapping;
global handles;
global num txt raw xls;
global overlapTransparent;
global h_grid_axis;

if monkey == 15
    M = 3;  % Recording + Microstimulation + Inactivation
elseif monkey == 13
    M =2;
elseif monkey == 10
    M= 1;     % Added for Messi @20220716
end

for m = 1:M   % Recording + Microstimulation + Inactivation
    
    hemisphere_text = {'Left','Right'};
    
    % Which channel do we select?
    
    pos = get(h_grid_axis,'CurrentPoint');
    pos = pos(1,1:2);
    xSelect = round(pos(2));
    ySelect = round(pos(1)- 0.5 * ~mod(xSelect,2));
    
    figure(802);
    
    delete(handles);
    
    % handles(100) = rectangle('Position',[gridRange(2,2)+1.4,xSelect-0.1,0.2,0.2],'FaceColor','k');
    handles(1) = plot(xlim,repmat(xSelect - 0.5 + overlapping(1,1),1,2),'r--','LineWid',1);
    handles(2) = plot(xlim,repmat(xSelect + 0.5 + overlapping(1,2),1,2),'r--','LineWid',1);
    handles(3) = plot(repmat(ySelect + 0.25,1,2),[xSelect - 0.5 + overlapping(1,1) xSelect + 0.5 + overlapping(1,2)],'r','LineWid',1);
    
    % Show corresponding sessions
    sessions_match = [];
    for channel = 1:length(data)
        if data{channel}{2}(1) == xSelect && data{channel}{2}(2) == ySelect
            sessions_match = [sessions_match data{channel}{1}];
        end
    end
    
    title(sprintf('Monkey %g, session(s) for %s [%g,%g]: %s',monkey,hemisphere_text{hemisphere},xSelect,ySelect,num2str(sessions_match)));
    
    
    figurePosition = [ 700+m*50  100  459  600 ];
    set(figure(802+m),'Position',figurePosition,'color','w'); clf
    
    
    % 2-D Visualization (Coronal)
    h_coronal = axes('Position',[0.15 0.1 0.8 0.8]);
    
    axis ij; hold on;
    
    % Overlapping MRI data
    try
        if hemisphere == 1 && monkey == 15
            fileNo = (xSelect - AP0) + Duncan_left_AP0 ;  % Because I use Polo_right MRI as standard images
        elseif monkey == 10
            fileNo = xSelect + Messi_left_AP0 - AP0;
        elseif monkey == 13
            fileNo = xSelect + Fara_left_AP0 - AP0;
            
            %     else
            %         fileNo = xSelect + Polo_right_AP0 - AP0;
        end
        
        MRI = imread([MRI_path num2str(fileNo) '.bmp']);
        
        % h_MRI = image(MRI_offset{1}+xSelect*MRI_offset{3}(1), MRI_offset{2} + xSelect*MRI_offset{3}(2),MRI);
        
        % MRI_offset_new = {[x_center, y_center, x_range], [slopes]}
        ratioyx = size(MRI,1)/size(MRI,2) * 10 * 0.8; % Auto keep ratio !!!   HH20180606
        x_lims = [MRI_offset_new{1}(1) - MRI_offset_new{1}(3)/2, MRI_offset_new{1}(1) + MRI_offset_new{1}(3)/2];
        y_lims = [MRI_offset_new{1}(2) - ratioyx * MRI_offset_new{1}(3)/2, MRI_offset_new{1}(2) + ratioyx * MRI_offset_new{1}(3)/2];
        
        if hemisphere==1
            MRI = flip(MRI,2);
        end
        h_MRI = image(x_lims + xSelect * MRI_offset_new{2}(1), y_lims + xSelect * MRI_offset_new{2}(2), MRI);
        
        set(h_MRI,'AlphaData',1);
        
        % 20180604
        uicontrol('style','pushbutton','unit','norm','pos', [0.0630    0.9510    0.0300    0.0210],...
            'callback',{@ manual_adjust_mri, 1});
        uicontrol('style','pushbutton','unit','norm','pos', [0.0630    0.9210    0.0300    0.0210],...
            'callback',{@ manual_adjust_mri, 2});
        uicontrol('style','pushbutton','unit','norm','pos', [0.0330    0.9350    0.0300    0.0210],...
            'callback',{@ manual_adjust_mri, 3});
        uicontrol('style','pushbutton','unit','norm','pos', [0.0930    0.9350    0.0300    0.0210],...
            'callback',{@ manual_adjust_mri, 4});
        
        uicontrol('style','pushbutton','unit','norm','pos', [0.1630    0.9510    0.0300    0.0210],...
            'callback',{@ manual_adjust_mri, 5});
        uicontrol('style','pushbutton','unit','norm','pos', [0.1630    0.9210    0.0300    0.0210],...
            'callback',{@ manual_adjust_mri, 6});
        
        uicontrol('style','pushbutton','unit','norm','pos', [0.2430    0.9510    0.0300    0.0210],...
            'callback',{@ manual_adjust_mri, 7});
        uicontrol('style','pushbutton','unit','norm','pos', [0.2430    0.9210    0.0300    0.0210],...
            'callback',{@ manual_adjust_mri, 8});
        
    catch
        disp('No MRI data found...');
    end
    
    % Frame
    xlim([gridRange(2,1) gridRange(2,2)+1]);
    ylim([-30 maxZ]);
    
    grid minor;
    set(h_coronal,'XMinorGrid','on','XMinorTick','on');
    % title(sprintf('Monkey %g, %s [%g, %g]',monkey, hemisphere,xSelect,ySelect));
    title(sprintf('Monkey %g, %s[%g], AP \\approx %g',monkey, hemisphere_text{hemisphere},xSelect,(AP0-xSelect)*0.8));
    
    % Keep scale
    aspectRatio = (range(ylim) * 100) / (range(xlim) * 800);  % grid interval = 0.8 mm
    set(figure(802+m),'Position',[figurePosition(1:2) figurePosition(4)/aspectRatio figurePosition(4)]);
    daspect([1/0.8 10 1]); % This is betther man
    
    for channel = 1:length(data)
        %     if data{channel}{2}(1) == xSelect  % Only plot the line we select
        if data{channel}{2}(1) >= xSelect + overlapping(1,1) && data{channel}{2}(1) <= xSelect + overlapping(1,2)   % Overlapping neighboring slices
            yLoc = data{channel}{2}(2)-0.5;
            GMData = data{channel}{4};
            if isempty(GMData); continue; end;
            
            GuideTubeAndOffset = data{channel}{3};
            offSet = round((GuideTubeAndOffset(2) + GuideTubeAndOffset(1) - 1.4) * 100);  % Related to Guide Tube 1.4 cm!!
            
            for GMType = -size(GMTypes,1):-1  % For each area type
                GMThisType = find(GMData(:,1) == GMType);   % Read out ranges for this type
                if isempty(GMThisType); continue; end       % If absent, next type
                
                for i = 1:length(GMThisType)  % For each appearance
                    zBegin = GMData(GMThisType(i),2) + offSet;
                    zEnd = GMData(GMThisType(i),3) + offSet;
                    
                    if overlapping(1,1)~=0 || overlapping(1,2)~=0  % If we overlap neighboring slices
                        p = patch([yLoc yLoc yLoc+1 yLoc+1],[zEnd,zBegin,zBegin,zEnd],GMTypes{-GMType,2},'EdgeColor','none','FaceAlpha',overlapTransparent);
                    else
                        rectangle('Position',[yLoc,zBegin,1,zEnd-zBegin],'LineWidth',linWid,...
                            'EdgeColor',GMTypes{-GMType,2});
                    end
                end
                
            end
            
            for GMType = -size(GMTypes,1):-1  % For each area type (that we are NOT SURE!!)
                GMThisType = find(GMData(:,1) == GMType - 100);   % Read out ranges for this type (that we are NOT SURE!!)
                if isempty(GMThisType); continue; end       % If absent, next type
                
                for i = 1:length(GMThisType)  % For each appearance
                    zBegin = GMData(GMThisType(i),2) + offSet;
                    zEnd = GMData(GMThisType(i),3) + offSet;
                    
                    % Note we use dotted line here to mark areas that we are not sure
                    rectangle('Position',[yLoc,zBegin,1,zEnd-zBegin],'LineWidth',linWid,'EdgeColor',GMTypes{-GMType,2},'LineStyle',':');
                end
                
            end
            
            % Add start and end markers
            if start_end_markers
                if yLoc+0.5 == ySelect
                    col_temp = 'c';
                    rectangle('Position',[yLoc,offSet,1,1],'LineWidth',linWid,'FaceColor',col_temp,'EdgeColor',col_temp);
                    offSet_selected = offSet;
                    if length(data{channel})<5
                        rectangle('Position',[yLoc,offSet + movDis,1,1],'LineWidth',linWid,'FaceColor',col_temp,'EdgeColor',col_temp);
                    else
                        rectangle('Position',[yLoc,offSet + data{channel}{5},1,1],'LineWidth',linWid,'FaceColor',col_temp,'EdgeColor',col_temp);
                    end
                else
                    rectangle('Position',[yLoc,offSet,1,1],'LineWidth',linWid,'FaceColor','k');
                    if length(data{channel})<5
                        rectangle('Position',[yLoc,offSet + movDis,1,1],'LineWidth',linWid,'FaceColor','k');
                    else
                        rectangle('Position',[yLoc,offSet + data{channel}{5},1,1],'LineWidth',linWid,'FaceColor','k');
                    end
                end
                
                
            end
            
        end
    end
    
    if exist('h_MRI')
        if exist('offSet_selected')
            set(h_MRI,'ButtonDownFcn',{@ShowDepth,offSet_selected});
        else
            set(h_MRI,'ButtonDownFcn','');
        end
    else
        set(gca,'color',[1 1 1]);
    end
    
    if start_end_markers
        % Channel indicator
        plot([ySelect ySelect],ylim,'r--','LineW',0.5);
    end
    
    xlabel('Grid Y No. (x 0.8 mm)');
    
    % Set y scale to mm
    ytick_temp = 0:50:maxZ;
    set(gca,'ytick',ytick_temp);
    set(gca,'yticklabel',ytick_temp/10);
    ylabel('Depth (mm)');
    
    % rectangle('Position',[ySelect+0.2,maxZ*0.95,0.6,10],'FaceColor','r','EdgeColor','r');
    
    if hemisphere == 1
        set(gca,'xdir','rev');
    end
    
    for oldsize = 5:100
        set(findall(gcf,'fontsize',oldsize),'fontsize',13);
    end
    set(findall(gcf,'tickdir','i'),'tickdir','o');
    
    drawnow;
    
end
%% Mark the recording sites
% figure(803);   % Overlaping in this figure
mark_color = [41 89 204; 248 28 83; 14 153 46]/255;
mark_shape = {'o', '^','*'};

%------------ CD SUs ----------------
% Mask: monkey & hemisphere & xLoc & significant choice preference
if monkey == 15
    load('Z:\Data\Tempo\Batch\CD_m15&13_SU_Heading\m15_selected_cell')
    monkey_selected_cell = m15_selected_cell';
    
    mask_tuning = {
        (strcmp(txt(:, xls.Protocol),'HD')) & (strcmpi(txt(:, xls.Area), 'CD')) ...
        & (num(:,xls.HD_rep) >= 8) & (num(:, xls.Units_RealSU) == 1) & (~strcmpi(txt(:, xls.Unit), 'MU21'))...
        & (num(:, xls.Monkey) == monkey) & num(:, xls.Hemisphere) == hemisphere ...
        & (num(:, xls.Xloc) >= xSelect + overlapping(2,1)) & (num(:, xls.Xloc) <= xSelect + overlapping(2,2));
        
        (strcmp(txt(:, xls.Protocol),'u-stim')) & (strcmpi(txt(:, xls.Area), 'CD')) ...
        & (num(:,xls.rep_ustim) >= 8) & (num(:, xls.uA) <= 80) & (strcmp(txt(:, xls.Note), '0:1500'))...
        & (num(:, xls.Monkey) == monkey) & num(:, xls.Hemisphere) == hemisphere ...
        & (num(:, xls.Xloc) >= xSelect + overlapping(2,1)) & (num(:, xls.Xloc) <= xSelect + overlapping(2,2));
        
        (strcmp(txt(:, xls.Protocol),'BHD')) & (strcmpi(txt(:, xls.Area), 'CD')) ...
        & (num(:,xls.Psy_n) >= 8) & (contains(txt(:, xls.Note), 'After-injection')) & (contains(txt(:, xls.Note), '0h'))...
        & (num(:, xls.Monkey) == monkey) & num(:, xls.Hemisphere) == hemisphere ...
        & (num(:, xls.Xloc) >= xSelect + overlapping(2,1)) & (num(:, xls.Xloc) <= xSelect + overlapping(2,2));
        };
    
elseif monkey ==13
    
    load('Z:\Data\Tempo\Batch\CD_m15&13_SU_Heading\m13_selected_cell')
    monkey_selected_cell = m13_selected_cell'; 
    
    mask_tuning = {
    (strcmp(txt(:, xls.Protocol),'HD')) & (strcmpi(txt(:, xls.Area), 'CD')) ...
    & (num(:,xls.HD_rep) >= 8) & (num(:, xls.Units_RealSU) == 1) & (~strcmpi(txt(:, xls.Unit), 'MU21'))...
    & (num(:, xls.Monkey) == monkey) & num(:, xls.Hemisphere) == hemisphere ...
    & (num(:, xls.Xloc) >= xSelect + overlapping(2,1)) & (num(:, xls.Xloc) <= xSelect + overlapping(2,2));   
    
    (strcmp(txt(:, xls.Protocol),'u-stim')) & (strcmpi(txt(:, xls.Area), 'CD')) ...
    & (num(:,xls.rep_ustim) >= 8) & (num(:, xls.uA) <= 80) & (strcmp(txt(:, xls.Note), '0:1500'))...
    & (num(:, xls.Monkey) == monkey) & num(:, xls.Hemisphere) == hemisphere ...
    & (num(:, xls.Xloc) >= xSelect + overlapping(2,1)) & (num(:, xls.Xloc) <= xSelect + overlapping(2,2));
    };

    
elseif monkey == 10
    mask_tuning = {
        (strcmp(txt(:, xls.Protocol),'HD')) & (strcmpi(txt(:, xls.Area), 'CD')) ...
        & (num(:,xls.HD_rep) >= 8) & (num(:, xls.Units_RealSU) == 1) & (~strcmpi(txt(:, xls.Unit), 'MU21'))...
        & (num(:, xls.Monkey) == monkey) & num(:, xls.Hemisphere) == hemisphere ...
        & (num(:, xls.Xloc) >= xSelect + overlapping(2,1)) & (num(:, xls.Xloc) <= xSelect + overlapping(2,2));
        };
end

for m = 1:size(mask_tuning,1)
    figure(802+m)
    to_plot = num(mask_tuning{m,1},:);
    
    if m == 1
        to_plot = to_plot(monkey_selected_cell==1,:); 
    end
    text(max(xlim)*0.6, max(ylim),sprintf('No. cells = %g',sum(mask_tuning{m,1})));
    
%     % Changed by ZZ to mark the cites according to the suggestions from GY
%     % 20230524
%     if m == 1    % Heading task 
%         choice_pref_all = to_plot(:,[xls.HD_vest_ChoicePref xls.HD_vis_ChoicePref xls.HD_comb_ChoicePref]);
%         % ========================== TODO ========================
%         % I need to change code in Cum_heading_discrim_ZZ to let ChoicePref
%         % has direction signals 
%         if_ipsilateral= to_plot(:, xls.Hemisphere) == hemisphere; 
%         choicepref_p_all = to_plot(:, [xls.HD_vest_ChoicePref_p xls.HD_vis_ChoicePref_p xls.HD_comb_ChoicePref_p]);
%         
%         [~,sort_ind]  = sort(choice_pref_all(:,1)); 
%         pos_ind = sign(choice_pref_all(:,1)) > 0; 
%         [~, pos_sort_ind] = sort(choice_pref_all(pos_ind,1));    % Let vestibular condition as example
%         [~, neg_sort_ind] = sort(choice_pref_all(~pos_ind, 1));
%         
% %         pos_colorcode = linspace(0,1, sum(pos_ind));            % (1,0,0) red for positve
% %         neg_colorcode = linspace(0,1, sum(~pos_ind));  % (0,0,1) blue for negtive
%         
%     elseif m == 2  % Microstimulation task 
%         dPSE_p_all = to_plot(:, [xls.PSE_p_vest xls.PSE_p_vis xls.PSE_p_comb]);
%         dPSE_all = to_plot(:, [xls.PSE_ctrl_vest xls.PSE_ctrl_vis xls.PSE_ctrl_comb]) - to_plot(:, [xls.PSE_stim_vest xls.PSE_stim_vis xls.PSE_stim_comb]);  % For Left hemisphere, positive value indicates contralateral shif
%         
%         [~,sort_ind]  = sort(dPSE_all(:,1)); 
%         pos_ind = sign(dPSE_all(:,1)) > 0; 
% %         [~, pos_sort_ind] = sort(dPSE_all(pos_ind,1));    % Let vestibular condition as example
% %         [~, neg_sort_ind] = sort(dPSE_all(~pos_ind, 1));
%         
% %         pos_colorcode = linspace(0,1, sum(pos_ind));            % (1,0,0) red for positve
% %         neg_colorcode = linspace(0,1, sum(~pos_ind));  % (0,0,1) blue for negtive
%     end
    
    for i = 1:size(to_plot,1)
        offSet = round((to_plot(i, xls.guidetube) + to_plot(i, xls.offset) - 1.4) * 100);
        %     offSet =  to_plot(i, xls.offset) ;
        xx = to_plot(i, xls.Yloc) + (rand(1) - 0.5)*0;
        yy = offSet + (to_plot(i, xls.Depth)/100);
        
        if m == 1
            % Plot according to the max choice divergence (if all non-significant, black circles)
            choice_pref_this = to_plot(i, xls.HD_vest_ChoicePref : xls.HD_comb_ChoicePref_p);
            [~, max_ind] = max(abs(choice_pref_this([1 3 5])));
            any_sig = any(choice_pref_this([2 4 6]) <= 0.05);
            if any_sig > 0
                mark_color_this = mark_color(max_ind,:); 
            else
                mark_color_this = [1 1 1];
            end
%         choice_pref_this = choice_pref_all(i, 1);
%         [~, max_ind] = max(abs(choice_pref_this))
%         any_sig = choicepref_p_all(i,1) < 0.05; 
%         if choice_pref_this > 0
%            where_is = find(choice_pref_all(pos_ind,1)==choice_pref_this);   % where this data is in positive dataset 
%            mark_color = [pos_colorcode(pos_sort_ind == where_is) 0 0];
%         else 
%            where_is = find(choice_pref_all(~pos_ind,1)==choice_pref_this);   % where this data is in negative dataset 
%            mark_color = [0 0 neg_colorcode((neg_sort_ind == where_is))];
%         end

        elseif m == 2
            dPSE_this = to_plot(i, [xls.PSE_ctrl_vest xls.PSE_ctrl_vis xls.PSE_ctrl_comb]) - to_plot(i, [xls.PSE_stim_vest xls.PSE_stim_vis xls.PSE_stim_comb]);
            % For Left hemisphere, positive value indicates contralateral shift
            [~,max_ind] = max(abs(dPSE_this));
            dPSE_p_all = to_plot(i, [xls.PSE_p_vest xls.PSE_p_vis xls.PSE_p_comb]);
            any_sig = any(dPSE_p_all <= 0.05);
            
            if any_sig > 0
                mark_color_this = mark_color(max_ind,:);
            else
                mark_color_this = [1 1 1];
            end
%             if dPSE_this > 0
%                 where_is = find(dPSE_all(pos_ind,1) == dPSE_this);
%                 mark_color = [pos_colorcode(pos_sort_ind == where_is) 0 0] ;
%             else
%                 where_is = find(dPSE_all(~pos_ind,1) == dPSE_this);
%                 mark_color = [0 0 neg_colorcode(flipud(neg_sort_ind) == where_is)] ;
%             end
            
%             if_contra = false;    % whether shift to contralateral side 
%             dPSE_p = to_plot(i, [xls.PSE_p_vest xls.PSE_p_vis xls.PSE_p_comb]);
%             any_sig = any(dPSE_p <= 0.05);
%             if any_sig
%                 if dPSE_p(2) <= 0.05
%                     max_ind = 2;
%                     dPSE = sign(to_plot(i,xls.PSE_ctrl_vis) - to_plot(i,xls.PSE_stim_vis)); 
%                 elseif dPSE_p(3) <= 0.05
%                     max_ind = 3;
%                     dPSE = sign(to_plot(i,xls.PSE_ctrl_comb) - to_plot(i,xls.PSE_stim_comb)); 
%                 else
%                     max_ind =1;
%                     dPSE = sign(to_plot(i,xls.PSE_ctrl_vest) - to_plot(i,xls.PSE_stim_vest));      % Postive indicate rightward shift 
%                 end
%                 if hemisphere == 1 && dPSE > -1
%                     if_contra = 1;
%                 elseif hemisphere == 2 && dPSE < 1
%                     if_contra =1;
%                 end
%             end
        elseif m ==3
            any_sig = false;
            mark_color_this = [0 0 0];
        end
        
%         if any_sig
%             if m ==1
                plot(xx,yy,mark_shape{m},'MarkerEdgeColor','k','linewid',0.1,'markerfacecol',mark_color_this,'markersize',7);
%                 plot(xx,yy,mark_shape{m},'MarkerEdgeColor','w','linewid',0.1,'markerfacecol',mark_color{max_ind},'markersize',7);
%             elseif m==2
%                 if if_contra
%                     plot(xx,yy,mark_shape{m},'MarkerEdgeColor',mark_color{max_ind},'linewid',0.1,'markerfacecol','none','markersize',7);
%                 else
%                     plot(xx,yy,mark_shape{m},'MarkerEdgeColor','k','linewid',0.1,'markerfacecol',mark_color_this,'markersize',7);
%                 end
%             end
%         else
%             plot(xx,yy,mark_shape{m},'MarkerEdgeColor',mark_color,'linewid',1,'markerfacecol','none','markersize',7);
%         end
        
        
    end
end

%% Plot Tuning.  HH20140624

%{
figure(803);
% --------------- Tuning Properties
% Mask: monkey & hemishpere & xLoc & significant visual tuning
% mask_tuning = (num(:,xls.Monkey) == monkey) & num(:,xls.Hemisphere)==hemisphere & (num(:,xls.Xloc) >= xSelect + overlapping(2,1) & num(:,xls.Xloc) <= xSelect + overlapping(2,2)) ...
%     & (num(:,xls.p_vis) < 0.05);
% to_plot = num(mask_tuning,:);
%
% for i =  1:size(to_plot,1)
%     offSet = round((to_plot(i,xls.guidetube)  + to_plot(i,xls.offset) - 2.0) * 100);  % Related to Guide Tube 2.0 cm!!
%     xx = to_plot(i,xls.Yloc) + ((to_plot(i,xls.Pref_vis) > 0) * 0.2 + (to_plot(i,xls.Pref_vis) <= 0) * -0.2)*sign((hemisphere==2)-0.5);  % Left and Right Prefer
%     yy = offSet + round(to_plot(i,xls.Depth)/100);
%
%     if to_plot(i,xls.Pref_vis) > 0
%         plot(xx,yy,'r>','linewid',1.2);
%     else
%         plot(xx,yy,'r<','linewid',1.2);
%     end
% %     plot(xx,yy,'ro' );
% end

% --------------- Memsac Properties
% % Mask: monkey & hemishpere & xLoc & significant visual tuning
% % mask_memsac = (num(:,xls.Monkey) == monkey) & (strcmp(txt(:,xls.Hemisphere),hemisphere)) & (num(:,xls.Xloc) == xSelect) ...
% %     & (num(:,xls.p_M) < 0.05);
% mask_memsac = (num(:,xls.Monkey) == monkey) & num(:,xls.Hemisphere)==hemisphere & (num(:,xls.Xloc) >= xSelect + overlapping(2,1) & num(:,xls.Xloc) <= xSelect + overlapping(2,2)) ...
%     & (1|num(:,xls.p_M) < 0.05) & (strcmp(txt(:,xls.Area),'LIP')) & (num(:,xls.Chan1)>0);
% to_plot = num(mask_memsac,:);
%
% for i =  1:size(to_plot,1)
%     offSet = round((to_plot(i,xls.guidetube)  + to_plot(i,xls.offset) - 2.0) * 100);  % Related to Guide Tube 2.0 cm!!
%     xx = to_plot(i,xls.Yloc) + ((to_plot(i,xls.pref_M) > 0) * 0.2 + (to_plot(i,xls.pref_M) <= 0) * -0.2)*sign((hemisphere==2)-0.5);  % Left and Right Prefer
%     yy = offSet + round(to_plot(i,xls.Depth)/100);
%     %
%     %     if to_plot(i,xls.pref_M) > 0
%     %         plot(xx,yy,'k>','linewid',1.2);
%     %     else
%     %         plot(xx,yy,'k<','linewid',1.2);
%     %     end
%     if to_plot(i,xls.p_M)<0.05
%         if to_plot(i,xls.pref_M) > 0
%             plot(xx,yy,'k>','linewid',1,'markerfacecol','g','markersize',6);
%         else
%             plot(xx,yy,'k<','linewid',1,'markerfacecol','g','markersize',6);
%             %            plot(xx,yy,'ko','linewid',1,'markerfacecol','k');
%         end
%
%     else
%         plot(to_plot(i,xls.Yloc),yy,'go','markersize',5);
%     end
% end

% --------------- MemSac Properties (MU & SU)
% Mask: monkey & hemishpere & xLoc & significant visual tuning
%{
mask_tuning = (num(:,xls.Monkey) == monkey) & num(:,xls.Hemisphere)==hemisphere & ...
    (num(:,xls.Xloc) >= xSelect + overlapping(2,1) & num(:,xls.Xloc) <= xSelect + overlapping(2,2)) ;
to_plot_num = num(mask_tuning,:);
to_plot_txt = txt(mask_tuning,:);

for i =  1:size(to_plot_num,1)
    
    if strcmp(to_plot_txt(i,xls.Protocol),'MemSac')
        offSet = round((to_plot_num(i,xls.guidetube)  + to_plot_num(i,xls.offset) - 2.0) * 100);  % Related to Guide Tube 2.0 cm!!
        xx = to_plot_num(i,xls.Yloc); % + ((to_plot(i,xls.Pref_vis) > 0) * 0.2 + (to_plot(i,xls.Pref_vis) <= 0) * -0.2)*sign((hemisphere==2)-0.5);  % Left and Right Prefer
        yy = offSet + round(to_plot_num(i,xls.Depth)/100);
        
        if to_plot_num(i,xls.p_M) < 0.01 %to_plot_num(i,xls.HD_MemSac) >= 0.8  % T site
            plot(xx,yy,'ro','linewid',0.4,'markerfacecol','r','markersize',5);
            pref_M = headingToAzi(to_plot_num(i,xls.pref_M))*pi/180;
            aspect = daspect;
            plot([xx xx+0.5*cos(pref_M)*sign((hemisphere==2)-0.5)],[yy yy+(-1)*0.5*sin(pref_M)*aspect(2)/aspect(1)],'r-','linew',1.5);
        else  % non-T site
            plot(xx,yy,'ro','linewid',0.4,'markerfacecol','none','markersize',5);
        end
    end
    
end
%}

% --------------- LIP SUs -------
% % Mask: monkey & hemishpere & xLoc & significant visual tuning
% %{

mask_tuning = (strcmp(txt(:,xls.Protocol),'HD') | strcmp(txt(:,xls.Protocol),'HD_dt')) & (strcmpi(txt(:,xls.Area),'LIP') | strcmpi(txt(:,xls.Area),'LIP-V')) ...
    & (num(:,xls.HD_rep) >= 8) & (num(:,xls.Units_RealSU) == 1) & num(:,xls.Monkey) == monkey & num(:,xls.Hemisphere)==hemisphere...
    & (num(:,xls.Xloc) >= xSelect + overlapping(2,1) & num(:,xls.Xloc) <= xSelect + overlapping(2,2)) & (num(:,xls.HD_TargFirst) == 1);

to_plot = num(mask_tuning,:);
text(max(xlim)*0.6, max(ylim),sprintf('No. cells = %g',sum(mask_tuning)));

for i =  1:size(to_plot,1)
    offSet = round((to_plot(i,xls.guidetube)  + to_plot(i,xls.offset) - 2.0) * 100);  % Related to Guide Tube 2.0 cm!!
    %    xx = to_plot(i,xls.Yloc) + ((to_plot(i,xls.Pref_vis) > 0) * 0.2 + (to_plot(i,xls.Pref_vis) <= 0) * -0.2)*sign((hemisphere==2)-0.5);  % Left and Right Prefer
    xx = to_plot(i,xls.Yloc) + (rand(1) - 0.5)*0;  % Left and Right Prefer
    yy = offSet + (to_plot(i,xls.Depth)/100);

    % Plot according to the max choice divergence (if all non-significant, black circles)
    max_color = {'b','r','g'};
    choice_pref_this = to_plot(i,xls.HD_vest_ChoicePref : xls.HD_comb_ChoicePref_p);
    [~, max_ind] = max(abs(choice_pref_this([1 3 5])));
    any_sig = any(choice_pref_this([2 4 6]) < 0.05);
    
    if any_sig
        plot(xx,yy,'wo','linewid',0.1,'markerfacecol',max_color{max_ind},'markersize',7);
    else
        plot(xx,yy,'ko','linewid',1,'markerfacecol','none','markersize',7);
%         plot(xx,yy,'ro','linewid',0.4,'markerfacecol','none','markersize',7);
    end

        
%     if to_plot(i,xls.HD_MemSac) >= 0.9 % T cell
%         plot(xx,yy,'ro','linewid',0.4,'markerfacecol','r','markersize',5);
%     else % non-T cell
%         plot(xx,yy,'ro','linewid',0.4,'markerfacecol','r','markersize',5);
% %         plot(xx,yy,'ro','linewid',0.4,'markerfacecol','none','markersize',7);
%     end
end
%}

function ShowDepth(~,~,offset)
persistent depthLine;
pos = get(gca,'CurrentPo');
depth = pos(1,2);

try
    set(depthLine(1),'ydata',[depth depth]);
    set(depthLine(2),'position',[max(xlim)*0.8 depth],'string',sprintf('%0.0f\n',(depth-offset)*100));
catch
    depthLine(1) = plot(xlim,[depth depth],'b--','ButtonDownFcn',{@ShowDepth,offset});
    depthLine(2) = text(max(xlim)*0.8,depth,sprintf('%0.0f\n',(depth-offset)*100),'color','b','fontsize',15);
end


% function ReadXls()
%
% %% Read xls for plotting tuning. HH20140624
% global num txt raw xls;
% [num,txt,raw] = xlsread('Z:\Data\MOOG\Results\Result.xlsm',2);
%
% % Get Header infomation
% HEADS_N = 3;
%
% header_all = txt(HEADS_N-1,:);
% header_all(strcmp(header_all,'')) = txt(HEADS_N-2,strcmp(header_all,''));
%
% for i = 1:length(header_all)
%     try
%         if i == num(1,i)
%             eval(['xls.' header_all{i} '=' num2str(i) ';']);
%         else
%             disp('Header info error...');
%             keyboard;
%         end
%     catch
%     end
% end
%
% % Delete headers
% end_line = find(~isnan(num(:,1)),1,'last');
% num = num(HEADS_N+1:end_line,:);
% txt = txt(HEADS_N : end_line - 1,:); % Note here
% raw = raw(HEADS_N+1:end_line,:);

function change_monkey_hemi(~,~)
DrawMapping_ZZ(get(gcbo,'value')); %gcbo：正在执行其回调的对象的句柄


% 20180604 Add manual align buttons
function manual_adjust_mri(source,event,choice)
global MRI_offset_new hemisphere;

manual_step_size = 0.5;
manual_zoom_size = 1.02;

switch choice
    case {1,2}
        MRI_offset_new{1}(2) = MRI_offset_new{1}(2) + sign(choice - 1.5) * manual_step_size * 10;
    case {3,4}
        MRI_offset_new{1}(1) = MRI_offset_new{1}(1) + sign(3.5 - choice) * sign(1.5 - hemisphere) * manual_step_size;
    case {5,6}
        MRI_offset_new{1}(3) = MRI_offset_new{1}(3) * manual_zoom_size^(sign(5.5 - choice));
    case {7,8}
        MRI_offset_new{2}(2) = MRI_offset_new{2}(2) + sign(choice - 7.5) * manual_step_size/3;
end

fprintf('MRI_offset_new = {[%g, %g, %g],[%g, %g]};\n',...
    MRI_offset_new{1}, MRI_offset_new{2});
SelectChannel(source,event);
