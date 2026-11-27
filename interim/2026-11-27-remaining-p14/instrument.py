"""Install the pre-trial width rule in the disposable p11 build only."""
import sys
from pathlib import Path

path = Path(sys.argv[1]) / "src/stages.rs"
src = path.read_text()
key = "        while progress < 1. - 1e-12 {\n            if !p8_parents.is_empty()"
assert src.count(key) == 1
src = src.replace(key, '''        let p14 = p8_reload && std::env::var_os("P14_GAP").is_some();
        let mut latest: Option<(f64, Vec<f64>)> = None;
        let mut prior: Option<(f64, Vec<f64>)> = None;
''' + key)
key = "            let next = if p8_predict || p9_window || local {"
assert src.count(key) == 1
src = src.replace(key, '''            let nominal = if p14 && progress >= 0.8401877295048414 {
                let pts = &committed.joints.as_ref().unwrap().points;
                let gaps: Vec<f64> = pts.iter().map(|p|
                    p.normal.iter().zip(p.jump_m).map(|(a,b)| a*b).sum()).collect();
                if latest.as_ref().is_none_or(|(t,_)| *t != progress) {
                    prior = latest.take();
                    latest = Some((progress,gaps.clone()));
                }
                let (left,right) = (0.9034196486885337,0.9034214389463122);
                let half = right-left;
                let mut proposal = nominal;
                let mut point = None;
                if progress < left && nominal >= left {
                    proposal = left;
                    point = Some(30);
                } else if left <= progress && progress < right {
                    proposal = nominal.min(progress+2.*half);
                    point = Some(30);
                }
                let mut forecasts = Vec::new();
                if let Some((t,old)) = &prior {
                    for (i,(&gap,&last)) in gaps.iter().zip(old).enumerate() {
                        if gap >= 0. || last >= gap { continue }
                        // Once point 30's measured window is reached, track the
                        // adjacent switches by their own committed secants.
                        if progress >= left && ![29,31].contains(&i) { continue }
                        let distance = -gap * (progress-t) / (gap-last);
                        let zero = progress + distance;
                        if zero <= nominal {
                            // A finite crossing width avoids iterating to zero
                            // from below. The ordinary acceptance audit still runs.
                            let to = progress + distance.max(2.*half);
                            forecasts.push((i,last,gap,zero));
                            if to < proposal { proposal = to; point = Some(i); }
                        }
                    }
                }
                eprintln!("{}",serde_json::json!({"p14_parent":progress,"attempt":count,
                    "previous":prior.as_ref().map(|(t,_)|t),"nominal":nominal,
                    "to":proposal,"point":point,"forecasts":forecasts}));
                proposal
            } else { nominal };
''' + key)
path.write_text(src)
