# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : extract_pdf_diagrams.R                                     ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 22-04-2025                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────
# Extract and crop diagrams from PDF pages
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────

library(magick)
library(pdftools)

#' Extract and crop figures from PDF pages
#'
#' @param input_path String: Path to the PDF file.
#' @param output_dir String: Directory to save PNGs.
#' @param pages Vector: Integer vector of pages to process (e.g., 2:15).
#' @param dpi Integer: Resolution (default 450 for high quality).
#' @param margin_top Integer: Pixels to cut from the top (removes titles).
#' @param margin_bottom Integer: Pixels to cut from the bottom (removes page numbers).
#' @param fuzz Integer: Sensitivity for smart trimming (0-100).
extract_pdf_figures <- function(input_path, 
                                output_dir, 
                                pages, 
                                dpi = 450, 
                                margin_top = 1000, 
                                margin_bottom = 1000,
                                fuzz = 15) {
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Iterate through pages
  for (p in pages) {
    message(sprintf("Processing page %d...", p))
    
    # 1. Read single page (memory efficient)
    img <- image_read_pdf(input_path, pages = p, density = dpi)
    info <- image_info(img)
    
    # 2. Calculate safe crop area
    crop_height <- info$height - margin_top - margin_bottom
    
    if (crop_height <= 0) {
      warning(sprintf("Page %d: Margins overlap. Skipping.", p))
      next
    }
    
    # 3. Define geometry: Keep full width, computed height, offset by top margin
    # Format: "WxH+X+Y"
    geo_string <- paste0(info$width, "x", crop_height, "+0+", margin_top)
    
    # 4. Processing pipeline: Crop Margins -> Smart Trim -> Save
    img_processed <- img %>%
      image_crop(geometry = geo_string) %>%
      image_trim(fuzz = fuzz)
    # image_transparent("white") # Uncomment if transparent background is needed
    
    # 5. Write to disk
    out_path <- file.path(output_dir, sprintf("figure_page_%02d.png", p))
    image_write(img_processed, path = out_path, format = "png")
  }
  
  message(sprintf("Done! Images saved in: %s", output_dir))
}

# Extract model diagrams from pdf file
extract_pdf_figures(
  input_path = "R scripts/Model diagrams/model_diagrams.pdf",
  output_dir = "Model diagrams",
  pages = 2:15,
  dpi = 450,
  margin_top = 1000,
  margin_bottom = 1000
)
