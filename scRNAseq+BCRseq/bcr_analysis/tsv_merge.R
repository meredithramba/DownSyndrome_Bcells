### To facilitate the amalgamation of various samples for subsequent analysis, 
### it is necessary to augment the resulting table with a column labeled sample_id 
### for sample identification

append_sample_id_column <- function(input_file, output_file, sample_id) {
  # Read the TSV file
  df <- read.delim(input_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  
  # Append the sample_id column
  df$sample_id <- sample_id
  
  # Write the updated dataframe to a new TSV file
  write.table(df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  message("Column 'sample_id' successfully added with value: ", sample_id)
}


append_sample_id_column(input_file = "DM001_db-pass.tsv",
                        output_file = "DM001_db-pass_ID.tsv",
                        sample_id = "DM001")

append_sample_id_column(input_file = "DM002_db-pass.tsv",
                        output_file = "DM002_db-pass_ID.tsv",
                        sample_id = "DM002")

append_sample_id_column(input_file = "DM003_db-pass.tsv",
                        output_file = "DM003_db-pass_ID.tsv",
                        sample_id = "DM003")

append_sample_id_column(input_file = "DM004_db-pass.tsv",
                        output_file = "DM004_db-pass_ID.tsv",
                        sample_id = "DM004")

append_sample_id_column(input_file = "DM005_db-pass.tsv",
                        output_file = "DM005_db-pass_ID.tsv",
                        sample_id = "DM005")

append_sample_id_column(input_file = "DM006_db-pass.tsv",
                        output_file = "DM006_db-pass_ID.tsv",
                        sample_id = "DM006")

append_sample_id_column(input_file = "DM007_db-pass.tsv",
                        output_file = "DM007_db-pass_ID.tsv",
                        sample_id = "DM007")

append_sample_id_column(input_file = "DM008_db-pass.tsv",
                        output_file = "DM008_db-pass_ID.tsv",
                        sample_id = "DM008")

append_sample_id_column(input_file = "DM009_db-pass.tsv",
                        output_file = "DM009_db-pass_ID.tsv",
                        sample_id = "DM009")

append_sample_id_column(input_file = "DM0011_db-pass.tsv",
                        output_file = "DM0011_db-pass_ID.tsv",
                        sample_id = "DM0011")

### Furthermore, to ensure uniqueness across samples, 
### sequence IDs and cell barcodes may be modified by appending the sample ID to them

add_text_to_column <- function(input_file, output_file, column_number, append_text) {
  # Read the TSV file
  df <- read.delim(input_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  
  # Ensure the row number is valid
  if (column_number < 1 || column_number > ncol(df)) {
    stop("Invalid row number.")
  }
  
  # Append text to each element in the specified row
  df[,column_number] <- paste0(df[,column_number], append_text)

  # Write back to a new TSV file
  write.table(df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  message("Text successfully added to column ", column_number)
}

#sequence ID
add_text_to_column(input_file = "DM001_db-pass_ID.tsv",
                   output_file = "DM001_db-pass_ID_labeled.tsv",
                   column_number = 1, 
                   append_text = "_DM001")

add_text_to_column(input_file = "DM002_db-pass_ID.tsv",
                   output_file = "DM002_db-pass_ID_labeled.tsv",
                   column_number = 1, 
                   append_text = "_DM002")

add_text_to_column(input_file = "DM003_db-pass_ID.tsv",
                   output_file = "DM003_db-pass_ID_labeled.tsv",
                   column_number = 1, 
                   append_text = "_DM003")

add_text_to_column(input_file = "DM004_db-pass_ID.tsv",
                   output_file = "DM004_db-pass_ID_labeled.tsv",
                   column_number = 1, 
                   append_text = "_DM004")

add_text_to_column(input_file = "DM005_db-pass_ID.tsv",
                   output_file = "DM005_db-pass_ID_labeled.tsv",
                   column_number = 1, 
                   append_text = "_DM005")

add_text_to_column(input_file = "DM006_db-pass_ID.tsv",
                   output_file = "DM006_db-pass_ID_labeled.tsv",
                   column_number = 1, 
                   append_text = "_DM006")

add_text_to_column(input_file = "DM007_db-pass_ID.tsv",
                   output_file = "DM007_db-pass_ID_labeled.tsv",
                   column_number = 1, 
                   append_text = "_DM007")

add_text_to_column(input_file = "DM008_db-pass_ID.tsv",
                   output_file = "DM008_db-pass_ID_labeled.tsv",
                   column_number = 1, 
                   append_text = "_DM008")

add_text_to_column(input_file = "DM009_db-pass_ID.tsv",
                   output_file = "DM009_db-pass_ID_labeled.tsv",
                   column_number = 1, 
                   append_text = "_DM009")

add_text_to_column(input_file = "DM0011_db-pass_ID.tsv",
                   output_file = "DM0011_db-pass_ID_labeled.tsv",
                   column_number = 1, 
                   append_text = "_DM0011")

##merge files
# Read the two TSV files
df1 <- read.delim("DM001_db-pass_ID_labeled.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
df2 <- read.delim("DM002_db-pass_ID_labeled.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
df3 <- read.delim("DM003_db-pass_ID_labeled.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
df4 <- read.delim("DM004_db-pass_ID_labeled.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
df5 <- read.delim("DM005_db-pass_ID_labeled.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
df6 <- read.delim("DM006_db-pass_ID_labeled.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
df7 <- read.delim("DM007_db-pass_ID_labeled.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
df8 <- read.delim("DM008_db-pass_ID_labeled.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
df9 <- read.delim("DM009_db-pass_ID_labeled.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
df11 <- read.delim("DM0011_db-pass_ID_labeled.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

# Combine them vertically
merged_df <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9, df11)

# Write back to a new TSV
write.table(merged_df, "BCR_DM001-DM011_data.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
                

