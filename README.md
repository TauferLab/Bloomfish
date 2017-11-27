# Overview
Bloomfish is a distributed-memory k-mer counting framework based on a single-node k-mer counting framework [Jellyfish](https://github.com/gmarcais/Jellyfish.git) and a MapReduce over MPI framework [Mimir](https://github.com/TauferLab/Mimir.git). It aims to combine advantages of both frameworks.

# Publication
* Tao Gao, Yanfei Guo, Yanjie Wei, Bingqiang Wang, Yutong Lu, Pietro Cicotti,
Pavan Balaji, and Michela Taufer. Bloomfish: A Highly Scalable Distributed K-mer
Counting Framework. IEEE International Conference on Parallel and Distributed
Systems (ICPADS) 2017.

# Installation
* git clone https://github.com/TauferLab/Bloomfish.git
* cd Bloomfish
* autoreconf -i
* ./configure --with-mimir=/mimir/install/dir
  --prefix=/bloomfish/install/dir
* make
* make install

# License

* The Mersenne Twister random generator is copyrighted by Agner Fog
  and distributed under the GPL version 3 or
  higher. http://www.agner.org.

* The Half float implementation is copyrighted by Industrial Light &
  Magic and is distributed under the license described in the
  HalfLICENSE file.

*   This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
