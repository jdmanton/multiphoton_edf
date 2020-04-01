%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% James Manton, 2020        %
% jmanton@mrc-lmb.cam.ac.uk %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SHOW_IMAGES = 1;
SAVE_IMAGES = 1;
ROOT_NAME = '0p65_water_5_steps_27um';

wavelength = 900E-9;           % Excitation wavelength (900 nm)
max_numerical_aperture = 0.70; % Numerical aperture of excitation lens
refractive_index = 1.33;       % Refractive index of focal volume
dof = 27E-6;                   % Target depth-of-field (27 microns)
num_steps = 5;                 % Number of zones in pupil mask
multiphoton_order = 2;         % 1 for single photon, 2 for two-photon, etc.


%% Calculate pupil zones to match DoF
nas = zeros([1, num_steps + 1]);
for i = 1:num_steps
    nas(i+1) = refractive_index * sind(acosd(cosd(asind(nas(i) / refractive_index)) - wavelength / (refractive_index * dof)));
end


%% Set up simulation parameters
M = 150;                % PSF planes will have dimensions (2M + 1) by (2M + 1)
N = 1315;               % Scaling factor -- increase for small voxel size
num_slices = 2 * M + 1; % Set so that PSF stack has same size in all dimensions
z_max = 43.9978E-06;    % Set so that Z range and hence voxel size is same as in XY

% Sampling params
m = -M:M;
n = m;
[m_2D, n_2D] = meshgrid(m, n);
k0 = 2 * pi / wavelength;
delta_k = k0 * max_numerical_aperture / M;
R = sqrt(m_2D.^2 + n_2D.^2);

% Aperture stop
stop = zeros(2 * M + 1, 2 * M + 1);
stop(R < M) = 1;

% Coordinates
theta = asin(delta_k  .* R / (k0 * refractive_index));
phi = atan2(n_2D, m_2D);
k_xy = delta_k * R;
k_z = sqrt((k0 * refractive_index)^2 - k_xy.^2);

% Distances
z_max = z_max / 2;
z_min = -z_max;
delta_x = M * wavelength / (max_numerical_aperture * N);
x_val = -delta_x * M : delta_x : delta_x * M;
y_val = x_val;
z_val = z_min : (z_max - z_min) / (num_slices - 1) : z_max;

disp(['refractive index: ', num2str(refractive_index)])
disp(['delta_x: ', num2str(delta_x)])
disp(['delta_z: ', num2str(z_val(2) - z_val(1))])
disp(['x FOV: ', num2str((x_val(end) - x_val(1)) * 1E6), ' microns'])
disp(['z FOV: ', num2str((z_max - z_min) * 1E6), ' microns'])


%% Pupil field
pupil_obliqueness = 1 ./ cos(theta);
pupil_polarisation(:, :, 1) = 1 / 2 * ((1 + cos(theta)) - (1 - cos(theta)) .* cos(2 * phi));
pupil_polarisation(:, :, 2) = -1 / 2 * (1 - cos(theta)) .* sin(2 * phi);
pupil_polarisation(:, :, 3) = -cos(phi) .* sin(theta);
pupil_phase = ones(size(pupil_obliqueness));

% Stop constraint
pupil_phase = pupil_phase .* stop;
pupil_polarisation = pupil_polarisation .* repmat(stop, [1, 1, 3]);
pupil_obliqueness = pupil_obliqueness .* stop;


%% Pupil amplitudes
r = R / max(m) * max_numerical_aperture;
intensity = zeros([size(pupil_obliqueness), num_slices, length(nas)]);
for i=1:(length(nas) - 1)
    pupil_amplitude = zeros(size(pupil_obliqueness));
    pupil_amplitude(r > nas(i) & r < nas(i+1)) = 1;
    pupil_amplitude = pupil_amplitude .* stop;

    % Focal field calculations
    electric_field = zeros([size(pupil_polarisation), num_slices]);
    progress = waitbar(0, ['Calculating for zone ', num2str(i), ' of ' num2str(num_steps), '...']);
    tic;
    for k_index = 1:num_slices
        electric_field(:, :, :, k_index) = debye.propagate(M, N, k_z, z_val(k_index), pupil_obliqueness, pupil_polarisation, pupil_amplitude, pupil_phase);
        intensity(:, :, k_index, i) = abs(electric_field(:, :, 1, k_index)).^2 + abs(electric_field(:, :, 2, k_index)).^2 + abs(electric_field(:, :, 3, k_index)).^2;
        intensity(:, :, k_index, i) = intensity(:, :, k_index, i).^multiphoton_order;
        waitbar(k_index / num_slices)
    end
    toc
    disp(['For ' num2str(num_slices) ' slices.'])
    close(progress)
end

%% Conventional focus
pupil_amplitude = zeros(size(pupil_obliqueness));
pupil_amplitude(r < nas(end)) = 1;
pupil_amplitude = pupil_amplitude .* stop;
electric_field = zeros([size(pupil_polarisation), num_slices]);
conventional_intensity = zeros([size(pupil_obliqueness), num_slices]);
progress = waitbar(0, 'Calculating conventional focus...');
tic;
for k_index = 1:num_slices
    electric_field(:, :, :, k_index) = debye.propagate(M, N, k_z, z_val(k_index), pupil_obliqueness, pupil_polarisation, pupil_amplitude, pupil_phase);
    conventional_intensity(:, :, k_index) = abs(electric_field(:, :, 1, k_index)).^2 + abs(electric_field(:, :, 2, k_index)).^2 + abs(electric_field(:, :, 3, k_index)).^2;
    conventional_intensity(:, :, k_index) = conventional_intensity(:, :, k_index).^multiphoton_order;
    waitbar(k_index / num_slices)
end
toc
disp(['For ' num2str(num_slices) ' slices.'])
close(progress)


%% Saving
if SAVE_IMAGES
	intensity_16bit = uint16(intensity / max(intensity(:)) * 2^16 - 1);
    for i = 1:(length(nas) - 1)
        for k_index = 1:num_slices
            imwrite(intensity_16bit(:, :, k_index, i), ['psfs/', ROOT_NAME, '-' num2str(i), '-', num2str(k_index, '%03d'), '-MP', num2str(multiphoton_order), '.tif'], 'tiff');
        end
    end
    conventional_intensity_16bit = uint16(conventional_intensity / max(conventional_intensity(:)) * 2^16 - 1);
    for k_index = 1:num_slices
        imwrite(conventional_intensity_16bit(:, :, k_index), ['psfs/', ROOT_NAME, '-' num2str(i), '-', num2str(k_index, '%03d'), '-MP', num2str(multiphoton_order), '-conv.tif'], 'tiff');
    end
end


%% Display
if SHOW_IMAGES
    % Zone PSFs
    for i=1:num_steps
        subplot(4, num_steps + 2, i)
        imshow(intensity(:, :, floor(size(intensity, 3) / 2), i), [0, max(intensity(:))])
        title(['Band ', num2str(i), ' XY'])
        
        subplot(4, num_steps + 2, num_steps + 2 + i)
        imshow(intensity(:, :, floor(size(intensity, 3) / 2), i), [0, max(intensity(:))])
        zoom(8)
        title(['Band ', num2str(i), ' XY zoomed'])
        
        subplot(4, num_steps + 2, 2 * (num_steps + 2) + i)
        imshow(squeeze(intensity(:, floor(size(intensity, 2) / 2), :, i))', [0, max(intensity(:))])
        title(['Band ', num2str(i), ' XZ'])
        
        subplot(4, num_steps + 2, 3 * (num_steps + 2) + i)
        imshow(squeeze(intensity(:, floor(size(intensity, 2) / 2), :, i))', [0, max(intensity(:))])
        zoom(8)
        title(['Band ', num2str(i), ' XZ zoomed'])
    end
    
    % Effective PSF
    subplot(4, num_steps + 2, num_steps + 1)
    imshow(sum(intensity(:, :, floor(size(intensity, 3) / 2), :), 4), [])
    title('Overall XY')
    
    subplot(4, num_steps + 2, 2 * (num_steps + 2) - 1)
    imshow(sum(intensity(:, :, floor(size(intensity, 3) / 2), :), 4), [])
    title('Overall XY')
    
    subplot(4, num_steps + 2, 3 * (num_steps + 2) - 1)
    imshow(sum(squeeze(intensity(:, floor(size(intensity, 2) / 2), :, :)), 3)', [])
    title('Overall XZ')
    
    subplot(4, num_steps + 2, 4 * (num_steps + 2) - 1)
    imshow(sum(squeeze(intensity(:, floor(size(intensity, 2) / 2), :, :)), 3)', [])
    zoom(8)
    title('Overall XZ zoomed')
    
    % Conventional PSF
    subplot(4, num_steps + 2, num_steps + 2)
    imshow(conventional_intensity(:, :, floor(size(conventional_intensity, 3) / 2)), [])
    title('Conventional XY')
    
    subplot(4, num_steps + 2, 2 * (num_steps + 2))
    imshow(conventional_intensity(:, :, floor(size(conventional_intensity, 3) / 2)), [])
    title('Conventional XY')
    
    subplot(4, num_steps + 2, 3 * (num_steps + 2))
    imshow(squeeze(conventional_intensity(:, floor(size(conventional_intensity, 2) / 2), :))', [])
    title('Conventional XZ')
    
    subplot(4, num_steps + 2, 4 * (num_steps + 2))
    imshow(squeeze(conventional_intensity(:, floor(size(conventional_intensity, 2) / 2), :))', [])
    zoom(8)
    title('Conventional XZ zoomed')
end
