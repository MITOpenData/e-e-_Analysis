./run.sh $1 "$2_beam_mix.root" $3 1 0 0 0 &>log/$2_beam_mix.log &
./run.sh $1 "$2_wta_mix.root" $3 1 0 1 0 &>log/$2_wta_mix.log &
./run.sh $1 "$2_thrust_mix.root" $3 1 1 0 0 &>log/$2_thrust_mix.log &

./run.sh $1 "$2_beam.root" $3 1 0 0 0 &>log/$2_beam.log &
./run.sh $1 "$2_wta.root" "" 1 0 1 0 &>log/$2_wta.log &
./run.sh $1 "$2_thrust.root" "" 1 1 0 0 &>log/$2_thrust.log &






