#' .DIANN_to_MQ_server2
#'
#' @description 
#' Server for app used to update PTM marks.
#' 
#' @param input The server function input, cf. https://shiny.posit.co/r/getstarted/shiny-basics
#' @param output The server function's output, cf. https://shiny.posit.co/r/getstarted/shiny-basics
#' @param session Shiny session.

.DIANN_to_MQ_server2 <- function(input,
                                 output,
                                 session) {
  homePath <- paste0(normalizePath(Sys.getenv("HOME"), winslash = "/"), "/R/proteoCraft")
  fl <- paste0(homePath, "/tmpMods.tmp") # We will use the same name always because we couldn't pass its name to the server.
  # Not being able to pass in-function (not global scope) variables to the server, despite it being instantiated within the function,
  # is the very reason I have to create this work-around.
  # Ultimately, this means we should not run parallel instance of the DIANN_to_MQ function.
  # But this is actually correct:
  #  - Who would want to resolve N+ shiny apps which opened in parallel!
  #  - The function takes a parallel cluster for internal calculations. We can't nest parallelization within parallelization.
  #
  con <- file(fl, "r+")
  ModTbl <- unserialize(connection = con)
  nr <- nrow(ModTbl)
  if (nr) {
    ModTbl$Name <- gsub(":", "_", ModTbl$Name) # This must be reverted before serializing!
    output$Msg <- shiny::renderUI({ shiny::em(" ") })
    output$PTMs <- shiny::renderUI({
      lst <- list()
      pick <- c()
      if (nr) {
        lst <- lapply(1L:nr, \(i) {
          pick[i] <- length(ModTbl$Options2[[i]]) >= 1L
          shiny::fluidRow(
            shiny::column(3L,
                          shiny::h5(shiny::strong(ModTbl$tag[i]))),
            shiny::column(4L,
                          shiny::textInput(paste0(ModTbl$Name[i], "_free"),
                                           NULL,
                                           ModTbl$Match[i])),
            shiny::column(4L,
                          if (pick[i]) {
                            #cat(ModTbl$Name[i], "\n")
                            shinyWidgets::pickerInput(ModTbl$Name[i],
                                                      NULL,
                                                      ModTbl$Options2[[i]],
                                                      ModTbl$Match[i],
                                                      FALSE,
                                                      shinyWidgets::pickerOptions(title = "Search me",
                                                                                  `live-search` = TRUE,
                                                                                  actionsBox = TRUE,
                                                                                  deselectAllText = "Clear search"))
                          } else {
                            shiny::em(paste0("(no match found in package unimod ", packageVersion("unimod"), ")"))
                          }
            )
            #,
            # shiny::column(1L,
            #               shiny::actionButton(paste0(ModTbl$Name[i], "_ok"), "ok"))
          )
        })
      }
      return(lst)
    })
    #
    sapply(1L:nr, \(i) {
      id <- ModTbl$Name[i]
      id2 <- paste0(id, "_free")
      if (length(ModTbl$Options2[[i]]) >= 1L) {
        shiny::observeEvent(input[[id]], {
          x <- as.character(input[[id]])
          #x <- UniMod2$Id[match(as.character(x), UniMod2$Id)] # In case of weird data conversion shenanigans
          if ((!is.na(x)) && nchar(x)) {
            shiny::updateTextInput(session,
                                   paste0(id, "_free"),
                                   NULL,
                                   x,
                                   x)
            ModTbl$Match[i] <- x
          }
          # if (sum(ModTbl$Match == "Nope...")) { shinyjs::disable("saveBtn") } else {
          #   shinyjs::enable("saveBtn")
          # }
        }, ignoreInit = TRUE)
      }
      shiny::observeEvent(input[[id2]], {
        x <- as.character(input[[id2]])
        if ((!is.na(x)) && nchar(x)) {
          ModTbl$Match[i] <<- x
        }
      })
      # shiny::observeEvent(input[[paste0(ModTbl$Name[i], "_ok")]], {
      #   id <- paste0(ModTbl$Name[i], "_free")
      #   #cat(input[[id]], "\n")
      #   if (input[[id]] != "") {
      #     ModTbl$Match[[i]] <- input[[id]]
      #   }
      #   if (sum(ModTbl$Match == "Nope...")) {  { shinyjs::disable("saveBtn") } else {
      #     shinyjs::enable("saveBtn")
      #   }
      # })
    })
    shiny::observeEvent(input$aaRemove, {
      sapply(1L:nr, \(i) {
        id2 <- paste0(ModTbl$Name[i], "_free")
        x <- gsub(":.*", "", as.character(input[[id2]]))
        if ((!is.na(x)) && nchar(x)) {
          shiny::updateTextInput(session,
                                 id2,
                                 NULL,
                                 x,
                                 x)
          ModTbl$Match[i] <<- x
        }
      })
    })
  }
  shiny::observeEvent(input$saveBtn, {
    suppressWarnings(suppressMessages(try({
      ModTbl$Name <- rownames(ModTbl)
      #Mods <<- ModTbl
      serialize(ModTbl, con)
      close(con)
    }, silent = TRUE))) # Only necessary in function mode
    shiny::stopApp()
  })
  # shiny::observeEvent(input$cancel, {
  #   suppressWarnings(suppressMessages(try({
  #     ModTbl$Name <- rownames(ModTbl)
  #     #Mods <<- ModTbl
  #     serialize(ModTbl, con)
  #     close(con)
  #   }, silent = TRUE))) # Only necessary in function mode
  #   shiny::stopApp()
  # })
  session$onSessionEnded(\() {
    suppressWarnings(suppressMessages(try({
      ModTbl$Name <- rownames(ModTbl)
      #Mods <<- ModTbl
      serialize(ModTbl, con)
      close(con)
    }, silent = TRUE))) # Only necessary in function mode
    shiny::stopApp()
  })
}
