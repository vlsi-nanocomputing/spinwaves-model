LIBRARY ieee;
USE ieee.std_logic_1164.all;
USE ieee.numeric_std.all;

ENTITY RCA_nbit IS
GENERIC(N : NATURAL:=32);
PORT
(
A : IN STD_LOGIC_VECTOR(N-1 DOWNTO 0);
B : IN STD_LOGIC_VECTOR(N-1 DOWNTO 0);
C : IN STD_LOGIC;
S : OUT STD_LOGIC_VECTOR(N DOWNTO 0)
);
END RCA_nbit;

ARCHITECTURE structural OF RCA_nbit IS

COMPONENT full_adder IS
PORT
(
A : IN STD_LOGIC;
B : IN STD_LOGIC;
carry_in : IN STD_LOGIC;
S : OUT STD_LOGIC;
C : OUT STD_LOGIC
);
END COMPONENT;

SIGNAL carry : STD_LOGIC_VECTOR(N DOWNTO 0);
------------------------ BEGIN ----------------------------------
BEGIN
carry(0) <= C;

GEN_FA:
FOR i IN 0 TO N-1 GENERATE
Fa: full_adder
PORT MAP(A(i), B(i), carry(i), S(i), carry(i+1));
END GENERATE GEN_FA;

S(N) <= carry(N);



END structural;
