-- Generated by Ora2Pg, the Oracle database Schema converter, version 17.4
-- Copyright 2000-2016 Gilles DAROLD. All rights reserved.
-- DATASOURCE: dbi:Oracle:host=sgd-nex2-db.stanford.edu;sid=SGD

SET client_encoding TO 'UTF8';

\set ON_ERROR_STOP ON


DROP TRIGGER IF EXISTS archcontigchange_bir ON nex.arch_contigchange CASCADE;
CREATE OR REPLACE FUNCTION trigger_fct_archcontigchange_bir() RETURNS trigger AS $BODY$
BEGIN
  IF (TG_OP = 'INSERT') THEN
 
       NEW.changed_by := upper(NEW.changed_by);
       PERFORM nex.checkuser(NEW.changed_by);

       RETURN NEW;
  END IF;

END;
$BODY$ LANGUAGE 'plpgsql';

CREATE TRIGGER archcontigchange_bir
BEFORE INSERT ON nex.arch_contigchange FOR EACH ROW
EXECUTE PROCEDURE trigger_fct_archcontigchange_bir();

DROP TRIGGER IF EXISTS archliteratureannotation_bir ON nex.arch_literatureannotation CASCADE;
CREATE OR REPLACE FUNCTION trigger_fct_archliteratureannotation_bir() RETURNS trigger AS $BODY$
BEGIN
  IF (TG_OP = 'INSERT') THEN

       NEW.created_by := upper(NEW.created_by);
       PERFORM nex.checkuser(NEW.created_by);

       RETURN NEW;
  END IF;

END;
$BODY$ LANGUAGE 'plpgsql';

CREATE TRIGGER archliteratureannotation_bir
BEFORE INSERT ON nex.arch_literatureannotation FOR EACH ROW
EXECUTE PROCEDURE trigger_fct_archliteratureannotation_bir();

DROP TRIGGER IF EXISTS archlocuschange_bir ON nex.arch_locuschange CASCADE;
CREATE OR REPLACE FUNCTION trigger_fct_archlocuschange_bir() RETURNS trigger AS $BODY$
BEGIN
  IF (TG_OP = 'INSERT') THEN

       NEW.added_by := upper(NEW.added_by);
       PERFORM nex.checkuser(NEW.added_by);

       RETURN NEW;
  END IF;

END;
$BODY$ LANGUAGE 'plpgsql';

CREATE TRIGGER archlocuschange_bir
BEFORE INSERT ON nex.arch_locuschange FOR EACH ROW
EXECUTE PROCEDURE trigger_fct_archlocuschange_bir();
