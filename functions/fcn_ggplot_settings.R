# ********************************************************
# GGPLOT settings ---------------------------------------
# ********************************************************

# Graphic specifications for ggplot:
themebar = theme(axis.text.x = element_text(face="bold", color="black", size=rel(1), angle=0),
                 axis.title.x = element_text(size = rel(1.1), angle = 0, face="bold"),
                 axis.text.y = element_text(face="bold", color="black", size=rel(1.05), angle=0), 
                 axis.title.y = element_text(size = rel(1.1), angle = 90, face="bold"),
                 plot.title = element_text(size = 16, angle = 0, face="bold"),
                 plot.subtitle = element_text(size = 14, angle = 0, face="bold"),
                 panel.border = element_rect(linetype = "solid", colour = "black", fill=NA),
                 legend.text = element_text(size = 12, face = "bold", lineheight=0.8, margin = margin(r = 0.3, unit = 'cm')),
                 legend.position = "bottom",
                 legend.box = "vertical",
                 legend.background = element_rect(fill=NA, linewidth=0.25, linetype="solid", colour ="black"),
                 legend.title = element_blank(),
                 legend.key = element_blank(),
                 legend.spacing.x=unit(0,"cm"), 
                 panel.grid.major = element_line(colour="gray", linetype = "dotted"),
                 panel.background = element_rect(fill = NA),
                 strip.background = element_rect(fill = NA),
                 strip.text = element_text(size=rel(1), face="bold"),
                 strip.placement = "outside")
