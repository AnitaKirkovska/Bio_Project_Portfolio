source("code/functions_edited.R")
folder <- "time3"

data <- extract_data_from_files(folder)
data <- trim_bad_spectra(data = data,bad_wavelengths = 180:229)
data <- remove_scans_with_bad_y_values(data = data,good_y_range = c(-200,100000))
data <- eliminate_overlap(data = data,dominant_spectra_method = "NIR",avg_by_id_first = TRUE)
# optional
data <- rescale_between_0_1(data,convert_to_absorbance = FALSE)
# smooth data to fill in missing data
data <- smooth_scans(data)

visualize_data(data)

# save data
# write_data_in_raw_read_format(data = data,filename = "")
# write_data_in_transposed_format(data = data,filename = "")