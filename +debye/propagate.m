function electric_field = propagate(M, N, kz, z, back_aperture_obliqueness, back_aperture_polarisation, back_aperture_amplitude, back_aperture_phase)

field_x = exp(1i * kz * z) .* back_aperture_obliqueness .* back_aperture_polarisation(:, :, 1) .* back_aperture_amplitude .* back_aperture_phase;
field_y = exp(1i * kz * z) .* back_aperture_obliqueness .* back_aperture_polarisation(:, :, 2) .* back_aperture_amplitude .* back_aperture_phase;
field_z = exp(1i * kz * z) .* back_aperture_obliqueness .* back_aperture_polarisation(:, :, 3) .* back_aperture_amplitude .* back_aperture_phase;

field_x = debye.paddedfft2(M, N, field_x);
field_y = debye.paddedfft2(M, N, field_y);
field_z = debye.paddedfft2(M, N, field_z);

electric_field = zeros([size(field_x), 3]);
electric_field(:, :, 1) = field_x;
electric_field(:, :, 2) = field_y;
electric_field(:, :, 3) = field_z;

end
