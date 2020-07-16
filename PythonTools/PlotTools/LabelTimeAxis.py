import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import datetime

def colour(dark):
    if dark:
        col = [0.7, 0.7, 0.7]
    else:
        col = [0.4, 0.4, 0.4]
    return col

def seconds(date):
    return (date - datetime.datetime.fromtimestamp(0)).total_seconds()

def LabelTimeAxis(ax, time, label_months = False, label_years = False):

    num_hour  = (time[-1] - time[0]) / (60*60)
    num_day   = num_hour / 24
    num_month = num_day / 30
    num_year  = num_day / 365

    draw_year  = 0
    draw_month = 0
    draw_day   = 0

    if num_hour < 500:
        draw_hour = 1
        if num_day > 1.5:
            draw_day = 1
    else:
        draw_hour = 0
        if num_day < 50:
            draw_day = 1
            if num_month > 1.5:
                draw_month = 1
        else:
            draw_day = 0
            if num_month < 50:
                draw_month = 1
                if num_year > 1.5:
                    draw_year = 1
            else:
                draw_month = 0
                draw_year = 1

    ax.set_ylim(0,draw_year+draw_month+draw_day+draw_hour)
    dark_year  = False
    dark_month = True
    dark_day   = False
    dark_hour  = True

    ticks = []
    labels = []
    if draw_year:
        ticks += [0.5]
        labels += ['Year']
    if draw_month:
        ticks += [draw_year+0.5]
        labels += ['Mo.']
    if draw_day:
        ticks += [draw_year+draw_month+0.5]
        labels += ['Day']
    if draw_hour:
        ticks += [draw_year+draw_month+draw_day+0.5]
        labels += ['Hour']
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels)

    one_day  = datetime.timedelta(0, 24 * 60 * 60)
    one_hour = datetime.timedelta(0,      60 * 60)

    # Indicate fast timescale
    start = datetime.datetime.fromtimestamp(time[ 0])
    stop  = datetime.datetime.fromtimestamp(time[-1]) 

    for year in range(start.year, stop.year+1):
    
        if draw_year == 1:
            left  = datetime.datetime(year,   1, 1)
            right = datetime.datetime(year+1, 1, 1)
    
            left  = seconds(max(left, start))
            right = seconds(min(right, stop))
    
            rect = mpl.patches.Rectangle( (left,0), right-left, 1, color=colour(dark_year))
            ax.add_patch(rect)
            dark_year = not(dark_year)
        
            if (label_years) or (draw_month == 1):
                ax.text( (left+right)/2, 0.5, str(year), ha='center', va='center')
    
        for month in range(1, 13):
        
            left  = datetime.datetime(year, month, 1, 0)
            if month < 12:
                right = datetime.datetime(year, month+1, 1, 0)
                num_days = (datetime.date(year, month+1, 1) - datetime.date(year, month, 1)).days
            else:
                right = datetime.datetime(year+1, 1, 1, 0)
                num_days = 31
        
            if (left < stop) and (right > start):
                
                if draw_month == 1:
                    left  = seconds(max(left, start))
                    right = seconds(min(right, stop))
    
                    rect = mpl.patches.Rectangle( (left, draw_year), right-left, 1, color=colour(dark_month))
                    ax.add_patch(rect)
                    dark_month = not(dark_month)
                
                    if (label_months) or (draw_year == 0):
                        month_str = datetime.date(1900, month, 1).strftime('%b')
                        ax.text( (left+right)/2, draw_year + 0.5, month_str, ha='center', va='center')
            
                for day in range(1, num_days+1):
                    left  = datetime.datetime(year, month, day, 0)
                    right = left + one_day
                
                    if (left < stop) and (right > start):
                    
                        if draw_day == 1:
                            left  = seconds(max(left, start))
                            right = seconds(min(right, stop))
                    
                            rect = mpl.patches.Rectangle( (left,draw_year+draw_month), right-left, 1, color=colour(dark_day))
                            ax.add_patch(rect)
                            dark_day = not(dark_day)
                            
                            if draw_month == 0:
                                ax.text( (left+right)/2, draw_year+draw_month+0.5, str(day), ha='center', va='center')
                
                        for hour in range(24):
                            left  = datetime.datetime(year, month, day, hour)
                            right = left + one_hour
                    
                            if (left < stop) and (right > start):
                
                                if draw_hour == 1:
                                    left  = seconds(max(left, start))
                                    right = seconds(min(right, stop))
        
                                    rect = mpl.patches.Rectangle( (left,draw_year+draw_month+draw_day), right-left, 1, color=colour(dark_hour))
                                    ax.add_patch(rect)
                                    dark_hour = not(dark_hour)


