#!/usr/bin/python
import argparse
import os
import numpy as np
import csv
import itertools
import math
import re

#GLOBALS
headers = ['time', 'frame', 'tunnel', 'fly_x', 'fly_y', 'fly_angle']
units = {
    'tunnel_dist' : 'mm',
    'chamber_dist' : 'mm',
    'tunnel_time' : 's',
    'chamber_time' : 's',
    'tunnel_vel' : 'mm/s',
    'chamber_vel' : 'mm/s',
    'tunnel_angular_vel' : 'degrees/s',
    'chamber_angular_vel' : 'degrees/s',
}

#CONSTANTS
CHAMBER_Y = 99
START_Y = 0

def find_experiments(directory):
    """
    Finds all flybar experiments for a directory.
    
    returns {
        date: [{name: "ExpName", path: "ExpPath"}, {name: "ExpName", path: "ExpPath"}],
        date: [{name: "ExpName", path: "ExpPath"}]
    }
    """
    datestamps = os.listdir(directory)
    exp_return = {}
    for d in datestamps:
        if not d.startswith("."):
            exp_return[d] = [] 
            full_d = os.path.join(directory,  d)
            exps = os.listdir(full_d)
            for e in exps:
                if not e.startswith("."):
                    exp_return[d].append({
                        'name' : e,
                        'path' : os.path.join(full_d, e),
                    })
    return exp_return

def load_data(exp, date):
    """
    Processes and analyzes data for a single experiment. A raw and analyzed output
    is written for the tunnel and chamber stages of the experiments.
    
    Returns 1 if loaded, 0 if did not load
    """
    trackers = []
    loaded = 0
    for fi in os.listdir(exp['path']):
        if fi.startswith("tracker"):
            trackers.append(fi)
        if not overwrite and fi.find("RAW") != -1:
            break
    #if loop is not broken
    else:
        run_dict = {}
        for run_file in sorted(trackers):
            p = re.compile('\d{6}')
            m = p.search(run_file)
            if not m:
                raise Error("Bad tracker name")
            run = m.group()
            tracker = load_tracking(os.path.join(exp['path'], run_file))
            run_dict[run] = parse_run(tracker)
        write_raw_output(run_dict, exp['path'], exp['name'], date)
        analyzed_dict = analyze_output(run_dict)
        write_analyzed_output(analyzed_dict, exp['path'], exp['name'], date)
        loaded = 1
    return loaded

def load_tracking(filepath):
    """
    Loads a single run into a numpy array, masks out bad data 
    returns a numpy array
    """
    #Cols: time,frame, tunnel, fly_x, fly_y, fly_angle, blob_ecc
    cols = [0,1,2,3,4,5,14]
    tracker_ar = np.loadtxt(filepath, delimiter="\t", usecols=cols, skiprows=1)
    mask =  (tracker_ar[:,6] != 0)
    return tracker_ar[mask]
    
def parse_run(fi_ar):
    """
    Takes the raw data from np array and splits out by tunnel and type
    output:
    return_dict = {
        tunnel 1 : {
            time: [...],
            frame: [...],
            fly_x: [...],
            fly_y: [...],
            fly_angle: [...]
        }
        tunnel 2 : { ... }
    }
    """
    header_ar = zip(range(0,len(headers)),headers)
    #removes the tunnel
    tunnel_row = header_ar[2]
    header_ar.remove(tunnel_row)
    lane_mask = []
    #mask out each lane (tunnel)
    for lane in range(0,6):
        lane_mask.append(fi_ar[:,2] == lane)
    return_dict = {}
    for tunnel, mask in enumerate(lane_mask):
        return_dict[tunnel+1] = {}
        for h in header_ar:
            return_dict[tunnel+1][h[1]] = fi_ar[:,h[0]][mask]
    return return_dict

def calculate_stats(x, y, time, angle):
    """
    For a single tunnel calculates the statistics for the tunnel and the etoh
    chamber.
    input: np arrays for x, y, time and angle for a single fly
    output: dictionary of statistics for tunnel and etoh chamber
    """
    if len(x) != len(y):
        pass
        #throw error
    stats_return = {
        'total_distance_tunnel' : 0,
        'total_distance_chamber' : 0,
        'total_time_tunnel' : 0,
        'total_time_chamber' : 0,
        'velocity_ar_tunnel' : [],
        'velocity_ar_chamber' : [],
        'angular_velocity_tunnel' : [],
        'angular_velocity_chamber' : [],
    }
    x0 = x.pop(0)
    y0 = y.pop(0)
    t0 = time.pop(0)
    a0 = angle.pop(0)
    for x1, y1, t1, a1 in itertools.izip(x, y, time, angle):
        if y1 <= START_Y:
            x0, y0, t0, a0 = x1, y1, t1, a1
            continue
        dist, t_delta, vel = calc_linear(x0,y0,t0,x1,y1,t1)
        angular_vel = calc_angular_velocity(a0, t0, a1, t1)
        if (y0 >= CHAMBER_Y and y1 >= CHAMBER_Y):
            stats_return['total_distance_chamber'] += dist
            stats_return['total_time_chamber'] += t_delta
            stats_return['velocity_ar_chamber'].append(vel)
            stats_return['angular_velocity_chamber'].append(angular_vel)
        elif ((START_Y < y0 <= CHAMBER_Y) and (START_Y < y1 <= CHAMBER_Y)): 
            stats_return['total_distance_tunnel'] += dist
            stats_return['total_time_tunnel'] += t_delta
            stats_return['velocity_ar_tunnel'].append(vel)
            stats_return['angular_velocity_tunnel'].append(angular_vel)
        x0, y0, t0, a0 = x1, y1, t1, a1
    return stats_return

def calc_angular_velocity(a0, t0, a1, t1):
    """
    Calculates the absolute value of angular velocity given two angles 
    and two time points. There is 180 degrees of ambiguity in a given angle.
    
    Angles are converted to radians and velocity is returned.
    """
    #detect if angles are closer forward or reverse
    delta_theta = math.fabs(a1-a0)
    if ((180 - delta_theta) < delta_theta):
        delta_theta = 180 - delta_theta
    angular_velocity = delta_theta/(t1-t0)
    return angular_velocity

def calc_linear(x0, y0, t0, x1, y1, t1):
    """
    Calculates linear statistics including change in time, distance covered
    and velocity.
    """
    t_delta = math.fabs(t1 - t0)
    dist = math.sqrt((x1 - x0)**2 + (y1 - y0) **2)
    vel = dist/t_delta
    return dist, t_delta, vel

def analyze_output(run_dict):
    """
    Takes the list of runs and tunnels, runs analysis on the values and
    combines them into a single data structure.
    """
    analyzed_dict = {}
    for run in run_dict:
        analyzed_dict[run] = {}
        analyzed_dict[run]['summary'] = {}
        for tunnel in run_dict[run]:
            analyzed_dict[run][tunnel] = {}
            x = run_dict[run][tunnel]['fly_x'].tolist()
            y = run_dict[run][tunnel]['fly_y'].tolist()
            time = run_dict[run][tunnel]['time'].tolist()
            angle = map(int, run_dict[run][tunnel]['fly_angle'].tolist())
            stats = calculate_stats(x, y, time, angle)
            #generate per tunnel dict
            analyzed_dict[run][tunnel]['tunnel_dist'] = stats['total_distance_tunnel']
            analyzed_dict[run][tunnel]['chamber_dist'] = stats['total_distance_chamber']
            analyzed_dict[run][tunnel]['tunnel_time'] = stats['total_time_tunnel']
            analyzed_dict[run][tunnel]['chamber_time'] = stats['total_time_chamber']
            analyzed_dict[run][tunnel]['tunnel_vel'] = np.mean(stats['velocity_ar_tunnel'])
            analyzed_dict[run][tunnel]['chamber_vel'] = np.mean(stats['velocity_ar_chamber'])
            analyzed_dict[run][tunnel]['tunnel_angular_vel'] = np.mean(stats['angular_velocity_tunnel'])
            analyzed_dict[run][tunnel]['chamber_angular_vel'] = np.mean(stats['angular_velocity_chamber'])
            if 'tunnel_dist' in analyzed_dict[run]['summary']:
                analyzed_dict[run]['summary']['tunnel_dist'] += stats['total_distance_tunnel']
                analyzed_dict[run]['summary']['chamber_dist'] += stats['total_distance_chamber']
                analyzed_dict[run]['summary']['tunnel_time'] += stats['total_time_tunnel']
                analyzed_dict[run]['summary']['chamber_time'] += stats['total_time_chamber']
                analyzed_dict[run]['summary']['tunnel_vel'].extend(stats['velocity_ar_tunnel'])
                analyzed_dict[run]['summary']['chamber_vel'].extend(stats['velocity_ar_chamber'])
                analyzed_dict[run]['summary']['tunnel_angular_vel'].extend(stats['angular_velocity_tunnel'])
                analyzed_dict[run]['summary']['chamber_angular_vel'].extend(stats['angular_velocity_chamber'])
            else:
                analyzed_dict[run]['summary']['tunnel_dist'] = stats['total_distance_tunnel']
                analyzed_dict[run]['summary']['chamber_dist'] = stats['total_distance_chamber']
                analyzed_dict[run]['summary']['tunnel_time'] = stats['total_time_tunnel']
                analyzed_dict[run]['summary']['chamber_time'] = stats['total_time_chamber']
                analyzed_dict[run]['summary']['tunnel_vel'] = stats['velocity_ar_tunnel']
                analyzed_dict[run]['summary']['chamber_vel'] = stats['velocity_ar_chamber']
                analyzed_dict[run]['summary']['tunnel_angular_vel'] = stats['angular_velocity_tunnel']
                analyzed_dict[run]['summary']['chamber_angular_vel'] = stats['angular_velocity_chamber']
        analyzed_dict[run]['summary']['tunnel_vel'] = np.mean(analyzed_dict[run]['summary']['tunnel_vel'])
        analyzed_dict[run]['summary']['chamber_vel'] = np.mean(analyzed_dict[run]['summary']['chamber_vel'])
        analyzed_dict[run]['summary']['tunnel_angular_vel'] = np.mean(analyzed_dict[run]['summary']['tunnel_angular_vel'])
        analyzed_dict[run]['summary']['chamber_angular_vel'] = np.mean(analyzed_dict[run]['summary']['chamber_angular_vel'])
    return analyzed_dict

def write_raw_output(run_dict, directory, experiment, date):
    """
    Takes the dictionary for all runs and condenses them to an output file.
    """
    tunnel_output = os.path.join(directory, "%s_%s_TUNNEL_RAW.txt" % (experiment,date))
    chamber_output = os.path.join(directory, "%s_%s_ENDCHAMBER_RAW.txt" % (experiment,date))
    tunnel_fi = open(tunnel_output, "w")
    chamber_fi = open(chamber_output, "w")
    output_headers = ['Experiment', 'Date', 'Time', 'Run', 'Tunnel', 'Parameter', 'Values']
    tunnel_fi.write("\t".join(output_headers))
    chamber_fi.write("\t".join(output_headers))
    tunnel_fi.write("\n")
    chamber_fi.write("\n")
    tcw = csv.writer(tunnel_fi, delimiter="\t")
    ccw = csv.writer(chamber_fi, delimiter="\t")
    for index, run in  enumerate(run_dict):
        for tunnel in run_dict[run]:
            mask_tunnel = (run_dict[run][tunnel]['fly_y'] < 99) & (run_dict[run][tunnel]['fly_y'] > 0)
            mask_chamber = run_dict[run][tunnel]['fly_y'] >= 99
            for header in headers:
                if header != 'tunnel':
                    tcw.writerow([experiment, date, run, index+1, tunnel, header] + run_dict[run][tunnel][header][mask_tunnel].tolist())
                    ccw.writerow([experiment, date, run, index+1, tunnel, header] + run_dict[run][tunnel][header][mask_chamber].tolist())
    tunnel_fi.close()
    chamber_fi.close()

def write_analyzed_output(analyzed_dict, directory, experiment, date):
    """
    Writes the analysis output file from the data structure.
    """
    tunnel_output = os.path.join(directory, "%s_%s_TUNNEL_ANALYZED.txt" % (experiment,date))
    chamber_output = os.path.join(directory, "%s_%s_ENDCHAMBER_ANALYZED.txt" % (experiment,date))
    tunnel_fi = open(tunnel_output, "w")
    chamber_fi = open(chamber_output, "w")
    analyzed_headers = ['Experiment', 'Date', 'Time', 'Run', 'Parameter', 'Units', 'Tunnel 1', 'Tunnel 2', 'Tunnel 3', 'Tunnel 4', 'Tunnel 5', 'Tunnel 6', 'Sum', 'Average']
    tunnel_fi.write("\t".join(analyzed_headers))
    chamber_fi.write("\t".join(analyzed_headers))
    tunnel_fi.write("\n")
    chamber_fi.write("\n")
    for index, run in enumerate(analyzed_dict):
        for parameter in sorted(analyzed_dict[run]['summary'].keys()):
            if parameter.find("chamber") != -1:
                ofi = chamber_fi
            else:
                ofi = tunnel_fi
            ofi.write("%s\t%s\t%s\t%s\t%s\t%s\t" % (experiment, date, run, index+1, parameter, units[parameter]))
            for tunnel in sorted(analyzed_dict[run].keys()):
                if parameter.find('vel') == -1 and tunnel == 'summary':
                    ofi.write("%0.5f\t%0.5f" % (analyzed_dict[run][tunnel][parameter], analyzed_dict[run][tunnel][parameter]/6.0))
                elif tunnel == 'summary':
                    ofi.write("%s\t%0.5f" % ('n/a', analyzed_dict[run][tunnel][parameter]))
                else:
                    ofi.write("%0.5f\t" % (analyzed_dict[run][tunnel][parameter]))
            ofi.write("\n")
    tunnel_fi.close()
    chamber_fi.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="Directory where tracking files are located")
    parser.add_argument("-v", "--verbose", dest='verbose', default=False, action="store_true", help="Flag: Increase output verbosity")
    parser.add_argument("-d", "--debug", dest='debug', default=False, action="store_true", help="Flag: Debugging option (also turns on verbosity)")
    parser.add_argument("--overwrite", dest='overwrite', default=False, action="store_true", help="Reanalyze all experiments and overwrite previous data (default is to only analyze new data)")
    args = parser.parse_args()
    
    #arguments
    directory = args.directory
    verbose = args.verbose
    debug = args.debug
    overwrite = args.overwrite
    if debug:
        verbose = True
    
    #find and load experiments
    experiments = find_experiments(directory)
    loaded = 0
    for date in experiments:
        for exp in experiments[date]:
            loaded += load_data(exp, date)
    
    #finished
    print "Experiments successfuly loaded: %d" % loaded

    
    
    
    
    