import subprocess

def merge_to_mp4(frame_filenames, movie_name, fps=12):
    f_log = open("ffmpeg.log", "w")
    f_err = open("ffmpeg.err", "w")
    cmd = ['ffmpeg', '-framerate', str(fps), '-i', frame_filenames, '-y', 
            '-q', '1', '-threads', '0', '-pix_fmt', 'yuv420p', movie_name]
    #print(cmd)
    subprocess.call(cmd, stdout=f_log, stderr=f_err)
    f_log.close()
    f_err.close()
