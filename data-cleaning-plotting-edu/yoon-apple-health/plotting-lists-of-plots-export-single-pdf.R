path = "~/Downloads/apple-health-dl/" 
zip = paste(path, 'export.zip', sep = '/')
unzip(zip, exdir = path)
Sys.sleep(3) # pause for 3 seconds to let your computer unzip it.
#list.files(paste0(path,'/apple_health_export'))

xml <- xmlParse(paste0(path, '/apple_health_export/export.xml'))
#summary(xml)
df_record <-   XML:::xmlAttrsToDataFrame(xml["//Record"]) %>% as_tibble()
df_activity <- XML:::xmlAttrsToDataFrame(xml["//ActivitySummary"]) %>% as_tibble()
df_workout <-  XML:::xmlAttrsToDataFrame(xml["//Workout"]) %>% as_tibble()

# df_clinical <- XML:::xmlAttrsToDataFrame(xml["//ClinicalRecord"])
# df_location <- XML:::xmlAttrsToDataFrame(xml["//Location"]) %>% 
#   mutate(latitude = as.numeric(as.character(latitude)),
#          longitude = as.numeric(as.character(longitude)))

df_record %>% names()
df = df_record %>%
  mutate(device = gsub(".*(name:)|,.*", "",device),
         value = as.numeric(as.character(value)),
         endDate = ymd_hms(endDate,tz="America/New_York"),
         date = date(endDate),
         year = year(endDate),
         month = month(endDate),
         day = day(endDate),
         yday = yday(endDate),
         wday = wday(endDate),
         hour = hour(endDate),
         minute = minute(endDate),
         type = str_remove(type, "HKQuantityTypeIdentifier")
  )

health_selection_list = df %>% select(type) %>% distinct()
  
# print(health_selection_list)
# health_selection_number = ask("What would you like to see? Select the number representing the item on the list: ")
# health_selection = health_selection_list[[1]][as.numeric(health_selection_number)]

# p1 = df %>%
#   filter(year(endDate) >= 2019) %>% 
#   arrange(endDate) %>% 
#   filter(type == health_selection) %>% 
#   # Had to reduce sourceName to these 2 sources to avoid double-counting
#   # filter(sourceName %in% c("Health")) %>% 
#   
#   ggplot(aes(x= date, y = value)) +
#   geom_point(alpha = 0.3) +
#   geom_smooth(span = 0.2, col = "grey30", se = FALSE) +
#   labs(title = health_selection
#     #  , caption = "200923"
#        ) 

#theme(axis.text.y = "Pounds")

plot_list = list()
  
for(i in 1:nrow(health_selection_list)){
  health_selection = health_selection_list[[1]][i]
  p1 = df %>%
    filter(year(endDate) >= 2019) %>% 
    arrange(endDate) %>% 
    filter(type == health_selection) %>% 
    # Had to reduce sourceName to these 2 sources to avoid double-counting
    # filter(sourceName %in% c("Health")) %>% 
    ggplot(aes(x= date, y = value)) +
    geom_point(alpha = 0.3) +
    geom_smooth(span = 0.2, col = "grey30", se = FALSE) +
    labs(title = health_selection
         ) 
  plot_list[[i]] = p1
  
  }
  
pdf("apple_health_custom_plots.pdf", onefile = TRUE)

for (i in seq(length(plot_list))) {
  print(plot_list[[i]])
}

dev.off()
