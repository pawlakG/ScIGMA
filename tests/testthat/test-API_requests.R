library(testthat)
library(ScIGMA)

test_that("fetch_clinical_vep_annotations handles successful response", {
  
  # Mock httr::POST to return a successful 200 response with some mock JSON
  mock_post <- function(url, ...) {
    response <- structure(list(), class = "response")
    response$status_code <- 200
    
    # Mock VEP JSON response
    mock_json <- '[
      {
        "input": "chr1:g.100>A",
        "transcript_consequences": [
          {
            "canonical": 1,
            "gene_symbol": "TEST_GENE",
            "hgvsp": "p.Ala123Val",
            "hgvsc": "c.123A>G",
            "transcript_id": "ENST000001",
            "variant_class": "SNV",
            "consequence_terms": ["missense_variant"]
          }
        ],
        "colocated_variants": []
      }
    ]'
    
    # Needs to match what httr::content returns for "parsed" simplifyVector=TRUE
    # But wait, httr::content parses the raw content using jsonlite.
    response$content <- charToRaw(mock_json)
    
    # We must also mock httr::content directly because parsing raw in tests can be tricky with content types
    return(response)
  }
  
  mock_content <- function(x, as, simplifyVector = FALSE, ...) {
    if (as == "text") return("Success")
    if (as == "parsed") {
      return(jsonlite::fromJSON('[
        {
          "input": "chr1:g.100>A",
          "transcript_consequences": [
            {
              "canonical": 1,
              "gene_symbol": "TEST_GENE",
              "hgvsp": "p.Ala123Val",
              "hgvsc": "c.123A>G",
              "transcript_id": "ENST000001",
              "variant_class": "SNV",
              "consequence_terms": ["missense_variant"]
            }
          ],
          "colocated_variants": []
        }
      ]', simplifyVector = simplifyVector))
    }
  }

  local_mocked_bindings(
    POST = mock_post,
    content = mock_content,
    .package = "httr"
  )
  
  # Call function (it maps chr1-100-G-A to chr1:g.100G>A usually, but our mock doesn't care about the exact POST body)
  # The function uses `original_variant = custom_variant_vector[match(input, hgvs_queries)]`
  # Our custom_variant_vector is "chr1-100-G-A", it translates to "chr1:g.100G>A"
  # Let's match the mock "input" to exactly what the function expects:
  mock_content_correct <- function(x, as, simplifyVector = FALSE, ...) {
    if (as == "text") return("Success")
    if (as == "parsed") {
      return(jsonlite::fromJSON('[
        {
          "input": "chr1:g.100G>A",
          "transcript_consequences": [
            {
              "canonical": 1,
              "gene_symbol": "TEST_GENE",
              "hgvsp": "p.Ala123Val",
              "hgvsc": "c.123A>G",
              "transcript_id": "ENST000001",
              "variant_class": "SNV",
              "consequence_terms": ["missense_variant"]
            }
          ]
        }
      ]', simplifyVector = simplifyVector))
    }
  }
  
  local_mocked_bindings(
    content = mock_content_correct,
    .package = "httr"
  )

  res <- ScIGMA:::fetch_clinical_vep_annotations(c("chr1-100-G-A"))
  
  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 1)
  expect_equal(res$gene, "TEST_GENE")
  expect_equal(res$protein, "TEST_GENE:p.A123V")
  expect_true(grepl("A123V", res$protein))
})

test_that("fetch_clinical_vep_annotations handles errors gracefully", {
  mock_post_error <- function(...) {
    response <- structure(list(), class = "response")
    response$status_code <- 500
    return(response)
  }
  
  mock_content_error <- function(x, as, ...) {
    return("Internal Server Error")
  }
  
  local_mocked_bindings(
    POST = mock_post_error,
    content = mock_content_error,
    .package = "httr"
  )
  
  expect_error(
    ScIGMA:::fetch_clinical_vep_annotations("chr1-100-G-A"),
    "API VEP 500: Internal Server Error"
  )
})
